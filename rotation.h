#ifndef ROTATION_H
#define ROTATION_H

#include <algorithm>
#include <cmath>
#include <limits>

//////////////////////////////////////////////////////////////////
// math functions needed for rotation conversion. 

// dot and cross production 
//内联函数定义点乘和叉乘，减小调用函数的开销
template<typename T>
inline T DotProduct(const T x[3], const T y[3]) {
    return (x[0] * y[0] + x[1] * y[1] + x[2] * y[2]);
}

template<typename T>
inline void CrossProduct(const T x[3], const T y[3], T result[3]) {
    result[0] = x[1] * y[2] - x[2] * y[1];
    result[1] = x[2] * y[0] - x[0] * y[2];
    result[2] = x[0] * y[1] - x[1] * y[0];
}


//////////////////////////////////////////////////////////////////


// 从旋转向量到四元素的变换theta*n->q=[q0,q1,q2,q3]
template<typename T>
inline void AngleAxisToQuaternion(const T *angle_axis, T *quaternion) {
    //起别名
    const T &a0 = angle_axis[0];
    const T &a1 = angle_axis[1];
    const T &a2 = angle_axis[2];
    //求theta，因为旋转向量是单位旋转轴乘上旋转角，所以先求出旋转角theta^2
    const T theta_squared = a0 * a0 + a1 * a1 + a2 * a2;

    //epsilon 是可以保证 1.0 + eps != 1.0 这个表达式成立的最小的正双精度浮点数。
    if (theta_squared > T(std::numeric_limits<double>::epsilon())) {
        const T theta = sqrt(theta_squared);
        //利用四元素到旋转向量的公式，q0=cos(theta/2),q1 = a0*sin(theta/2),...
        const T half_theta = theta * T(0.5);
        const T k = sin(half_theta) / theta;
        quaternion[0] = cos(half_theta);
        quaternion[1] = a0 * k;
        quaternion[2] = a1 * k;
        quaternion[3] = a2 * k;
    } else { 
        //这里theta为零，k = sin(theta/2)/theta 泰勒展开就是1/2
        const T k(0.5);
        quaternion[0] = T(1.0);
        quaternion[1] = a0 * k;
        quaternion[2] = a1 * k;
        quaternion[3] = a2 * k;
    }
}

//从四元数到旋转向量
template<typename T>
inline void QuaternionToAngleAxis(const T *quaternion, T *angle_axis) {
    const T &q1 = quaternion[1];
    const T &q2 = quaternion[2];
    const T &q3 = quaternion[3];
    const T sin_squared_theta = q1 * q1 + q2 * q2 + q3 * q3;

    //同样判断sin(theta)^2,这里的theta是旋转向量theta/2
    if (sin_squared_theta > T(std::numeric_limits<double>::epsilon())) {
        const T sin_theta = sqrt(sin_squared_theta);
        //得到cos(theta)，这里的theta同样是旋转向量theta的一半
        const T &cos_theta = quaternion[0];

        // If cos_theta is negative, theta is greater than pi/2, which
        // means that angle for the angle_axis vector which is 2 * theta
        // would be greater than pi...
        
        //由于cos(theta)为负数,而atan2(y,x)是四象限计算，范围在【-180，180】
        //那计算出来就会是>pi/2,导致轴角theta>pi所以要限制角度
        const T two_theta = T(2.0) * ((cos_theta < 0.0)
                                      ? atan2(-sin_theta, -cos_theta)
                                      : atan2(sin_theta, cos_theta));
        const T k = two_theta / sin_theta;

        angle_axis[0] = q1 * k;
        angle_axis[1] = q2 * k;
        angle_axis[2] = q3 * k;
    } else {
        // For zero rotation, sqrt() will produce NaN in derivative since
        // the argument is zero. By approximating with a Taylor series,
        // and truncating at one term, the value and first derivatives will be
        // computed correctly when Jets are used..
        const T k(2.0);//这里同样是通过；泰勒展开，这里是k=2theta/sin(theta) = 2
        angle_axis[0] = q1 * k;
        angle_axis[1] = q2 * k;
        angle_axis[2] = q3 * k;
    }

}

//求解旋转后的点，result = R*point,point在函数中使用pt表示
//R = cos(theta)I +(1-cos(theta))n(n)T,sin(theta)n^
//则result = cos(theta)pt+sin(theta)nxpt+(1-cos(theta))n(n)Tpt，(n)Tpt,实际就是n点乘p
template<typename T>
inline void AngleAxisRotatePoint(const T angle_axis[3], const T pt[3], T result[3]) {
    //求(theta)^2 = a0^2+a1^2+a2^2
    const T theta2 = DotProduct(angle_axis, angle_axis);
    if (theta2 > T(std::numeric_limits<double>::epsilon())) {
        // Away from zero, use the rodriguez formula
        //
        //   result = pt costheta +
        //            (w x pt) * sintheta +
        //            w (w . pt) (1 - costheta)
        //
        // We want to be careful to only evaluate the square root if the
        // norm of the angle_axis vector is greater than zero. Otherwise
        // we get a division by zero.
        //
        const T theta = sqrt(theta2);
        const T costheta = cos(theta);
        const T sintheta = sin(theta);
        const T theta_inverse = 1.0 / theta;
        //a0,a1,a2单位化，n0,n1,n2
        const T w[3] = {angle_axis[0] * theta_inverse,
                        angle_axis[1] * theta_inverse,
                        angle_axis[2] * theta_inverse};

        // Explicitly inlined evaluation of the cross product for
        // performance reasons.
        /*const T w_cross_pt[3] = { w[1] * pt[2] - w[2] * pt[1],
                                  w[2] * pt[0] - w[0] * pt[2],
                                  w[0] * pt[1] - w[1] * pt[0] };*/
        //n叉乘pt
        T w_cross_pt[3];
        CrossProduct(w, pt, w_cross_pt);

        //tmp为nTpt*(1-cos(theta)),nT为1x3,pt为3x1,所以tep是一个常量
        const T tmp = DotProduct(w, pt) * (T(1.0) - costheta);
        //    (w[0] * pt[0] + w[1] * pt[1] + w[2] * pt[2]) * (T(1.0) - costheta);
        //根据最开始写的公式带入即可
        result[0] = pt[0] * costheta + w_cross_pt[0] * sintheta + w[0] * tmp;
        result[1] = pt[1] * costheta + w_cross_pt[1] * sintheta + w[1] * tmp;
        result[2] = pt[2] * costheta + w_cross_pt[2] * sintheta + w[2] * tmp;
    } else {
        //对于零的情况这里解释的很清楚
        // Near zero, the first order Taylor approximation of the rotation
        // matrix R corresponding to a vector w and angle w is
        //
        //   R = I + hat(w) * sin(theta)
        //
        // But sintheta ~ theta and theta * w = angle_axis, which gives us
        //
        //  R = I + hat(w)
        //
        // and actually performing multiplication with the point pt, gives us
        // R * pt = pt + w x pt.
        //
        // Switching to the Taylor expansion near zero provides meaningful
        // derivatives when evaluated using Jets.
        //
        // Explicitly inlined evaluation of the cross product for
        // performance reasons.
        /*const T w_cross_pt[3] = { angle_axis[1] * pt[2] - angle_axis[2] * pt[1],
                                  angle_axis[2] * pt[0] - angle_axis[0] * pt[2],
                                  angle_axis[0] * pt[1] - angle_axis[1] * pt[0] };*/
        T w_cross_pt[3];
        CrossProduct(angle_axis, pt, w_cross_pt);

        result[0] = pt[0] + w_cross_pt[0];
        result[1] = pt[1] + w_cross_pt[1];
        result[2] = pt[2] + w_cross_pt[2];
    }
}

#endif // rotation.h