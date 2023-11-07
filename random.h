#ifndef RAND_H
#define RAND_H

#include <math.h>
#include <stdlib.h>

inline double RandDouble()
{   
    //强制转换,rand()的范围是从【0，RAND——MAX】
    double r = static_cast<double>(rand());
    return r / RAND_MAX;
}

/*利用均匀分布构造正态分布N(0,1)随机点
Marsaglia polar method算法：  得到满足N（0,1）正态分布的随机点
  和y是圆： x*x+y*y <=1内的随机点
  RandNormal()处理后  x*w是服从标准正态分布的样本
*/
inline double RandNormal()
{
    double x1, x2, w;
    do{
        x1 = 2.0 * RandDouble() - 1.0;
        x2 = 2.0 * RandDouble() - 1.0;
        w = x1 * x1 + x2 * x2;
    }while( w >= 1.0 || w == 0.0);

    w = sqrt((-2.0 * log(w))/w);
    return x1 * w;
}

#endif // random.h