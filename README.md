## Abstract

Camera calibration for estimating the intrinsic parameters and lens distortion is a prerequisite for various monocular vision applications including feature tracking and video stabilization. This application paper proposes a model for estimating the parameters on the fly by fusing gyroscope and camera data, both readily available in modern day smartphones. The model is based on joint estimation of visual feature positions, camera parameters, and the camera pose, the movement of which is assumed to follow the movement predicted by the gyroscope. Our model assumes the camera movement to be free, but continuous and differentiable, and individual features are assumed to stay stationary. The estimation is performed online using an extended Kalman filter, and it is shown to outperform existing methods in robustness and insensitivity to initialization. We demonstrate the method using simulated data and empirical data from an iPad.

## Paper

[arXiv](https://arxiv.org/abs/1805.12506)

## Video

<iframe width="560" height="315" src="https://www.youtube.com/embed/ro7TeQKgfT0" frameborder="0" gesture="media" allow="encrypted-media" allowfullscreen></iframe>

## Codes

[Codes on GitHub](https://github.com/AaltoVision/camera-gyro-calibration)

## Referencing

{% raw  %}
```
@inproceedings{Cortes+Solin+Kannala:2018,
      title = {Robust Gyroscope-Aided Camera Self-Calibration},
     author = {Cort{\'e}s, Santiago and Solin, Arno
               and Kannala, Juho},
       year = {2018},
  booktitle = {Proceedings of the International Conference on 
               Information Fusion (FUSION)}
}
```
{% endraw  %}
