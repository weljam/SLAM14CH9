#ifndef SnavelyReprojection_H
#define SnavelyReprojection_H

#include <iostream>
#include "ceres/ceres.h"
#include "rotation.h"

class SnavelyReprojectionError {
public:
    //初始化
    SnavelyReprojectionError(double observation_x, double observation_y) : observed_x(observation_x),
                                                                           observed_y(observation_y) {}
    //括号重载，仿函数
    template<typename T>
    bool operator()(const T *const camera,
                    const T *const point,
                    T *residuals) const {
        // camera是旋转角
        T predictions[2];
        CamProjectionWithDistortion(camera, point, predictions);
        //误差=投影像素点-观测像素点
        residuals[0] = predictions[0] - T(observed_x);
        residuals[1] = predictions[1] - T(observed_y);

        return true;
    }

    //camera：9维向量
    //[0-2]：角度轴旋转
    //[3-5]：平移
    //[6-8]：相机参数，[6]焦距，[7-8]是二阶和四阶径向畸变
    //point：地图点。
    //predictions：像素点。
    template<typename T>
    //静态函数编译时不会加入this指针，所以对于所有对象都共享
    static inline bool CamProjectionWithDistortion(const T *camera, const T *point, T *predictions) {
        //利用轴角求解point旋转后的点
        T p[3];
        AngleAxisRotatePoint(camera, point, p);
        //p再加上平移量得到变换后的点坐标
        p[0] += camera[3];
        p[1] += camera[4];
        p[2] += camera[5];

        //归一化平面，使用的数据是假设在投影平面之后，所以在去掉深度信息后，还需要乘以-1
        // Compute the center fo distortion
        T xp = -p[0] / p[2];
        T yp = -p[1] / p[2];

        //引用二阶和四阶径向畸变系数
        const T &l1 = camera[7];
        const T &l2 = camera[8];

        //利用切向畸变模型 xdis = x(1+k1r^2+k2r^4) 
        T r2 = xp * xp + yp * yp;
        T distortion = T(1.0) + r2 * (l1 + l2 * r2);
        
        //focal焦点
        const T &focal = camera[6];
        //u = fx*x_dis+cx,vfy*y_dis+cy
        predictions[0] = focal * distortion * xp;
        predictions[1] = focal * distortion * yp;

        return true;
    }

    static ceres::CostFunction *Create(const double observed_x, const double observed_y) {
        //输出维度2为误差项，即像素误差，输入维度，9为camera，3为地图点
        return (new ceres::AutoDiffCostFunction<SnavelyReprojectionError, 2, 9, 3>(
            new SnavelyReprojectionError(observed_x, observed_y)));
    }

private:
    double observed_x;
    double observed_y;
};

#endif // SnavelyReprojection.h

