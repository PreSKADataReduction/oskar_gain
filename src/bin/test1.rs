use hdf5::{File, Result};

use num::{
    complex::Complex
};
use ndarray::{
    ArrayView1
    , Array3
    
};

use lds::{
    station::{
        Station
    }
    , cfg::StationCfg
};

use serde_yaml::{
    from_reader
};


fn main() -> Result<()>{
    let station_cfg:StationCfg=from_reader(std::fs::File::open("../lds/rand_lfaa.yaml").unwrap()).unwrap();
    let station=Station::<Complex<f64>, f64>::from_cfg(&station_cfg);
    
    let ff:Vec<_>=station.fine_ch_freq_in_fs().into_iter().map(|f|{
        f/station.dt
    }).collect();
    let fc:Vec<_>=station.coarse_freq_of_fine_ch_in_fs().into_iter().map(|f|{
        f/station.dt
    }).collect();
    
    //let nants=station.ants.len();
    let nants=29;
    let nfreq=fc.len();
    let ntime=1;
    assert_eq!(fc.len(), ff.len());
    println!("{} {} {}", ff[0], ff[1]-ff[0], nfreq);
    let file = File::create("gain_model.h5")?; // open for writing
    let group=file.group("/")?;
    group.new_dataset_builder().with_data(ArrayView1::from(&ff)).create("freq (Hz)")?;

    let mut gain=unsafe{Array3::<[f64;2]>::uninit((ntime,nfreq,nants)).assume_init()};
    gain.iter_mut().for_each(|x|{
        x[0]=1.0;
        x[1]=0.0;
    });
    group.new_dataset_builder().with_data(gain.view()).create("gain_xpol")?;
    group.new_dataset_builder().with_data(gain.view()).create("gain_ypol")?;

    Ok(())
}

