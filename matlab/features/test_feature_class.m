clear all
close all
datapath_calib = '../data/desk';
[data, frames]=read_camera_data(datapath_calib,0);

feat=features(frames,100,[]);
feat.track_all();