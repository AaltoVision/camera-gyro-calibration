# camera-gyro-calibration
Santiago Cortés Reina (Doctoral student at Aalto University)


Summary
--
Camera autocalibration using the gyroscope. Using a linearized dynamic model and a pinhole projection camera model with two distortion coefficients.


Dependencies
--
Created in MatLab R2017a.

FFMPEG is used to extract the frames.

In order to extract features mexopencv (3.0.0) is used.

The extracted points are available and thus the example can be run  without FFMPEG and mexopencv.

Code description
--
All files are written in Mathworks Matlab

* `main.m` (Matlab)
  Example of selfcalibration
* `data`(folder)
  Folder containing example video, IMU data and extracted features.


Example 
--

* run `main.m` to generate the example solution as shown in 1.
* visualization control variables are present( video frames need to be extracted into data folder beforehand)
* to extract frames run the following comment.

```
ffmpeg -r 2 -i data/cards2/senserve-ios-2017-09-15-184221.mov data/cards2/frame-%05d.png 
```


License
--

This software is distributed under the GNU General Public License (version 3 or later); please refer to the file `LICENSE.txt`, included with the software, for details. 
