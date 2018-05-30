# camera-gyro-calibration
Santiago Cortés Reina (Doctoral student at Aalto University)


Summary
--
Camera autocalibration using the gyroscope. Using a linearized dynamic model and a pinhole projection camera model with two distortion coefficients.


Dependencies
--
FFMPEG.
In order to extract features mexopencv is used. However, the extracted points are available and thus it can be run 
as is.

Code description
--
All files are written in Mathworks Matlab

* `main.m` (Matlab)
  Example of selfcalibration
* `data`(folder)
  Folder containing example video, IMU data and extracted features.


How to generate the solution
--

* run `main.m` to generate the example solution.
* visualization control variables are present(caution, video frames are extracted into data folder)



License
--

This software is distributed under the GNU General Public License (version 3 or later); please refer to the file `LICENSE.txt`, included with the software, for details. 
