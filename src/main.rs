use hdf5::{File as H5File, Result, 
    H5Type,
};

use num::{
    complex::Complex
    , traits::FloatConst
};

use std::{
    fs::{
        create_dir_all
        , File
    }
    , io::Write
    , path::{
        PathBuf
    }
};

use clap::{
    Command
    , Arg
};

use ndarray::{
    ArrayView1
    , Array3
    
};

use lds::{
    station::{
        Station
    }
    , constants::light_speed
    , cfg::StationCfg
};

use serde_yaml::{
    from_reader
};

type FloatType=f32;
//const FloatSize:hdf5::types::FloatSize=hdf5::types::FloatSize::U8;

#[derive(H5Type, Clone, PartialEq, Debug)] // register with HDF5
#[repr(C)]
struct MyComplex{
    pub re: FloatType,
    pub im: FloatType
}

fn main() -> Result<()>{
    let matches=Command::new("oskar_station")
    .arg(
        Arg::new("cfg")
        .short('c')
        .long("cfg")
        .takes_value(true)
        .value_name("cfg file")
        .required(true)
    )
    .arg(
        Arg::new("out_dir")
        .short('o')
        .long("out")
        .takes_value(true)
        .value_name("output directory name")
        .required(true)
    )
    .arg(
        Arg::new("az0")
        .short('A')
        .long("az0")
        .takes_value(true)
        .value_name("pointing az0 in deg")
        .required(true)
    )
    .arg(
        Arg::new("ze0")
        .short('Z')
        .long("ze0")
        .takes_value(true)
        .value_name("pointing ze0 in deg")
        .required(true)
    )
    .get_matches();

    let station_cfg:StationCfg=from_reader(std::fs::File::open(matches.value_of("cfg").unwrap()).unwrap()).unwrap();
    let station=Station::<Complex<FloatType>, FloatType>::from_cfg(&station_cfg);

    let out_dir=PathBuf::from(matches.value_of("out_dir").unwrap());
    create_dir_all(&out_dir).unwrap();

    let mut layout_file=File::create(out_dir.join("layout.txt")).unwrap();
    for a in &station.ants{
        writeln!(&mut layout_file, "{},{},{}", a.pos[0], a.pos[1], a.pos[2]).unwrap();
    }

    
    let ff_Hz:Vec<_>=station.fine_ch_freq_in_fs().into_iter().map(|f|{
        f/station.dt
    }).collect();
    let fc_Hz:Vec<_>=station.coarse_freq_of_fine_ch_in_fs().into_iter().map(|f|{
        f/station.dt
    }).collect();

    let fc_fs=station.coarse_freq_of_fine_ch_in_fs();

    let az0=matches.value_of("az0").unwrap().parse::<FloatType>().unwrap().to_radians();
    let ze0=matches.value_of("ze0").unwrap().parse::<FloatType>().unwrap().to_radians();
    

    let delay=station.calc_required_digital_delay(az0, ze0);

    println!("{:?}", delay);
    
    
    //let nants=station.ants.len();
    let nants=station.ants.len();
    let nfreq=fc_Hz.len();
    let ntime=1;
    assert_eq!(fc_Hz.len(), ff_Hz.len());
    println!("{} {} {}", ff_Hz[0], ff_Hz[1]-ff_Hz[0], nfreq);
    let file = H5File::create(out_dir.join("gain_model.h5"))?; // open for writing
    let group=file.group("/")?;


    group.new_dataset_builder().with_data(ArrayView1::from(&ff_Hz)).create("freq (Hz)")?;

    let mut gain=unsafe{Array3::<MyComplex>::uninit((ntime,nfreq,nants)).assume_init()};
    
    for it in 0..ntime{
        for ifreq in 0..nfreq{
            for ia in 0..nants{
                let d=delay[ia];
                let g=Complex::<FloatType>::new(0.0, -2.0 * FloatType::PI() * fc_fs[ifreq] * d).exp();
                gain[(it, ifreq, ia)]=MyComplex{
                    re:g.re,im:g.im
                };
                //println!("{}", (2.0 * FloatType::PI() * fc_fs[ifreq] * d).to_degrees());
            }
        }
    }

    group.new_dataset_builder().with_data(gain.view()).create("gain_xpol")?;
    group.new_dataset_builder().with_data(gain.view()).create("gain_ypol")?;

    Ok(())
}

