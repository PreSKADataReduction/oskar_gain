#![allow(non_snake_case)]
use hdf5::{File as H5File, H5Type, Result};

use num::{complex::Complex, traits::FloatConst};

use std::{
    fs::{create_dir_all, File},
    io::Write,
    path::PathBuf,
};

use clap::Parser;

use ndarray::{Array3, ArrayView1};

use lds::{cfg::StationCfg, station::Station};

use serde_yaml::from_reader;

type FloatType = f32;
//const FloatSize:hdf5::types::FloatSize=hdf5::types::FloatSize::U8;

#[derive(Parser, Debug)]
#[clap(author, version, about)]
struct Args {
    #[clap(short('c'), long("cfg"), value_name("cfg file"))]
    cfg: String,

    #[clap(short('o'), long("out"), value_name("output directory name"))]
    outdir: String,

    #[clap(short('A'), long("az0"), value_name("pointing az0 in deg"))]
    azimuth0: FloatType,

    #[clap(short('Z'), long("ze0"), value_name("poiting ze0 in deg"))]
    zenith0: FloatType,
}

#[derive(H5Type, Clone, PartialEq, Debug)] // register with HDF5
#[repr(C)]
struct MyComplex {
    pub re: FloatType,
    pub im: FloatType,
}

fn main() -> Result<()> {
    let args = Args::parse();

    let station_cfg: StationCfg = from_reader(std::fs::File::open(args.cfg).unwrap()).unwrap();
    let station = Station::<Complex<FloatType>, FloatType>::from_cfg(&station_cfg);

    let out_dir = PathBuf::from(args.outdir);
    create_dir_all(&out_dir).unwrap();

    let mut layout_file = File::create(out_dir.join("layout.txt")).unwrap();
    for a in &station.ants {
        writeln!(&mut layout_file, "{},{},{}", a.pos[0], a.pos[1], a.pos[2]).unwrap();
    }

    let ff_Hz: Vec<_> = station
        .fine_ch_freq_in_fs()
        .into_iter()
        .map(|f| f / station.dt)
        .collect();
    let fc_Hz: Vec<_> = station
        .coarse_freq_of_fine_ch_in_fs()
        .into_iter()
        .map(|f| f / station.dt)
        .collect();

    let fc_fs = station.coarse_freq_of_fine_ch_in_fs();

    let az0 = args.azimuth0.to_radians();
    let ze0 = args.zenith0.to_radians();

    let delay = station.calc_required_digital_delay(az0, ze0);

    println!("{:?}", delay);

    //let nants=station.ants.len();
    let nants = station.ants.len();
    let nfreq = fc_Hz.len();
    let ntime = 1;
    assert_eq!(fc_Hz.len(), ff_Hz.len());
    println!("{} {} {}", ff_Hz[0], ff_Hz[1] - ff_Hz[0], nfreq);
    let file = H5File::create(out_dir.join("gain_model.h5"))?; // open for writing
    let group = file.group("/")?;

    group
        .new_dataset_builder()
        .with_data(ArrayView1::from(&ff_Hz))
        .create("freq (Hz)")?;

    let mut gain = unsafe { Array3::<MyComplex>::uninit((ntime, nfreq, nants)).assume_init() };

    for it in 0..ntime {
        for ifreq in 0..nfreq {
            for ia in 0..nants {
                let d = delay[ia];
                let g =
                    Complex::<FloatType>::new(0.0, -2.0 * FloatType::PI() * fc_fs[ifreq] * d).exp();
                gain[(it, ifreq, ia)] = MyComplex { re: g.re, im: g.im };
                //println!("{}", (2.0 * FloatType::PI() * fc_fs[ifreq] * d).to_degrees());
            }
        }
    }

    group
        .new_dataset_builder()
        .with_data(gain.view())
        .create("gain_xpol")?;
    group
        .new_dataset_builder()
        .with_data(gain.view())
        .create("gain_ypol")?;

    Ok(())
}
