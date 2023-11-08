#include <cstdio>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <Eigen/Core>
#include <Eigen/Dense>

#include "common.h"
#include "rotation.h"
#include "random.h"

/*
Map定义

Map(PointerArgType dataPtr, Index rows, Index cols, const StrideType& stride = StrideType())
可以看出,构建map变量,需要三个信息:指向数据的指针,构造矩阵的行数和列数

map相当于引用普通的c++数组,进行矩阵操作,而不用copy数据
*/
typedef Eigen::Map<Eigen::VectorXd> VectorRef;
typedef Eigen::Map<const Eigen::VectorXd> ConstVectorRef;

//fptr为文件流，format为读取格式，往value中添加数据
template<typename T>
void FscanfOrDie(FILE *fptr, const char *format, T *value) {
    int num_scanned = fscanf(fptr, format, value);
    if (num_scanned != 1)
        std::cerr << "Invalid UW data file. ";
}
//给point添加正态分布随机噪声
void PerturbPoint3(const double sigma, double *point) {
    for (int i = 0; i < 3; ++i)
        point[i] += RandNormal() * sigma;
}

//获得排在data中间的点，nth_element即排序第K小的值
double Median(std::vector<double> *data) {
    int n = data->size();
    std::vector<double>::iterator mid_point = data->begin() + n / 2;
    std::nth_element(data->begin(), mid_point, data->end());
    return *mid_point;
}

BALProblem::BALProblem(const std::string &filename, bool use_quaternions) {
    FILE *fptr = fopen(filename.c_str(), "r");

    if (fptr == NULL) {
        std::cerr << "Error: unable to open file " << filename;
        return;
    };

    //例：BAL数据集说明。第一行：
    //16 22106 83718
    //16个相机，22106个点，共进行83718次相机对点的观测
    //第2行到83719行：
    //6 18595     3.775000e+01 4.703003e+01
    //第6个相机观测18595个点，得到的相机的观测数据为3.775000e+01 4.703003e+01
    //第83720行到83720+16*9=83864
    //共16个相机的9纬参数：-R（3维），t（3维），f（焦距），k1,k2畸变参数
    //第83864到83864+3*22106=150182
    //22106个点的三维坐标

    FscanfOrDie(fptr, "%d", &num_cameras_);     //读取相机数量
    FscanfOrDie(fptr, "%d", &num_points_);      //读取路标点
    FscanfOrDie(fptr, "%d", &num_observations_);//读取观测点

    std::cout << "Header: " << num_cameras_
              << " " << num_points_
              << " " << num_observations_;

    //这里每一个观测点很多，所以以观测点的个数作为数组的大小，每一个观测点对应一个数组和一个相机
    point_index_ = new int[num_observations_];         //指针指向创建的路标点数组
    camera_index_ = new int[num_observations_];        //指针指向创建的相机数组
    observations_ = new double[2 * num_observations_]; //指针指向观测点数组

    num_parameters_ = 9 * num_cameras_ + 3 * num_points_;  //得到相机(9数据)+路标点(3数据)的总内存
    parameters_ = new double[num_parameters_];              //指针指向创建的总参数数组

    for (int i = 0; i < num_observations_; ++i) {
        FscanfOrDie(fptr, "%d", camera_index_ + i); //第几个相机
        FscanfOrDie(fptr, "%d", point_index_ + i);  //第几个路标点
        for (int j = 0; j < 2; ++j) {
            FscanfOrDie(fptr, "%lf", observations_ + 2 * i + j); //观测点
        }
    }

    for (int i = 0; i < num_parameters_; ++i) {
        FscanfOrDie(fptr, "%lf", parameters_ + i);  //读取相机的九个参数，轴角，平移，焦点，畸变参数
    }

    fclose(fptr);
    //是否使用四元数
    use_quaternions_ = use_quaternions;
    if (use_quaternions) {
        // 讲轴角变换为四元数
        num_parameters_ = 10 * num_cameras_ + 3 * num_points_;
        //扩大数组的大小
        double *quaternion_parameters = new double[num_parameters_];
        //相机参数和路标点的总参数
        double *original_cursor = parameters_;
        double *quaternion_cursor = quaternion_parameters;
        for (int i = 0; i < num_cameras_; ++i) {
            //在总参数中取出轴角，计算四元数，轴角+3，四元数+4
            AngleAxisToQuaternion(original_cursor, quaternion_cursor);
            quaternion_cursor += 4;
            original_cursor += 3;
            for (int j = 4; j < 10; ++j) {
                *quaternion_cursor++ = *original_cursor++;
            }
        }
        // Copy the rest of the points.
        for (int i = 0; i < 3 * num_points_; ++i) {
            *quaternion_cursor++ = *original_cursor++;
        }
        // Swap in the quaternion parameters.
        delete[]parameters_;
        parameters_ = quaternion_parameters;
    }
}

void BALProblem::WriteToFile(const std::string &filename) const {
    FILE *fptr = fopen(filename.c_str(), "w");

    if (fptr == NULL) {
        std::cerr << "Error: unable to open file " << filename;
        return;
    }

    fprintf(fptr, "%d %d %d %d\n", num_cameras_, num_cameras_, num_points_, num_observations_);

    for (int i = 0; i < num_observations_; ++i) {
        fprintf(fptr, "%d %d", camera_index_[i], point_index_[i]);
        for (int j = 0; j < 2; ++j) {
            fprintf(fptr, " %g", observations_[2 * i + j]);
        }
        fprintf(fptr, "\n");
    }

    for (int i = 0; i < num_cameras(); ++i) {
        double angleaxis[9];
        if (use_quaternions_) {
            //OutPut in angle-axis format.
            QuaternionToAngleAxis(parameters_ + 10 * i, angleaxis);
            memcpy(angleaxis + 3, parameters_ + 10 * i + 4, 6 * sizeof(double));
        } else {
            memcpy(angleaxis, parameters_ + 9 * i, 9 * sizeof(double));
        }
        for (int j = 0; j < 9; ++j) {
            fprintf(fptr, "%.16g\n", angleaxis[j]);
        }
    }

    const double *points = parameters_ + camera_block_size() * num_cameras_;
    for (int i = 0; i < num_points(); ++i) {
        const double *point = points + i * point_block_size();
        for (int j = 0; j < point_block_size(); ++j) {
            fprintf(fptr, "%.16g\n", point[j]);
        }
    }

    fclose(fptr);
}

// Write the problem to a PLY file for inspection in Meshlab or CloudCompare
void BALProblem::WriteToPLYFile(const std::string &filename) const {
    std::ofstream of(filename.c_str());

    of << "ply"
       << '\n' << "format ascii 1.0"
       << '\n' << "element vertex " << num_cameras_ + num_points_
       << '\n' << "property float x"
       << '\n' << "property float y"
       << '\n' << "property float z"
       << '\n' << "property uchar red"
       << '\n' << "property uchar green"
       << '\n' << "property uchar blue"
       << '\n' << "end_header" << std::endl;

    // Export extrinsic data (i.e. camera centers) as green points.
    double angle_axis[3];
    double center[3];
    for (int i = 0; i < num_cameras(); ++i) {
        const double *camera = cameras() + camera_block_size() * i;
        CameraToAngelAxisAndCenter(camera, angle_axis, center);
        of << center[0] << ' ' << center[1] << ' ' << center[2]
           << " 0 255 0" << '\n';
    }

    // Export the structure (i.e. 3D Points) as white points.
    const double *points = parameters_ + camera_block_size() * num_cameras_;
    for (int i = 0; i < num_points(); ++i) {
        const double *point = points + i * point_block_size();
        for (int j = 0; j < point_block_size(); ++j) {
            of << point[j] << ' ';
        }
        of << " 255 255 255\n";
    }
    of.close();
}

void BALProblem::CameraToAngelAxisAndCenter(const double *camera,
                                            double *angle_axis,
                                            double *center) const {
    VectorRef angle_axis_ref(angle_axis, 3);
    if (use_quaternions_) {
        QuaternionToAngleAxis(camera, angle_axis);
    } else {
        angle_axis_ref = ConstVectorRef(camera, 3);
    }

    // c = -R't
    Eigen::VectorXd inverse_rotation = -angle_axis_ref;
    AngleAxisRotatePoint(inverse_rotation.data(),
                         camera + camera_block_size() - 6,
                         center);
    VectorRef(center, 3) *= -1.0;
}

void BALProblem::AngleAxisAndCenterToCamera(const double *angle_axis,
                                            const double *center,
                                            double *camera) const {
    ConstVectorRef angle_axis_ref(angle_axis, 3);
    if (use_quaternions_) {
        AngleAxisToQuaternion(angle_axis, camera);
    } else {
        VectorRef(camera, 3) = angle_axis_ref;
    }

    // t = -R * c
    AngleAxisRotatePoint(angle_axis, center, camera + camera_block_size() - 6);
    VectorRef(camera + camera_block_size() - 6, 3) *= -1.0;
}

void BALProblem::Normalize() {
    // Compute the marginal median of the geometry
    std::vector<double> tmp(num_points_);
    Eigen::Vector3d median;
    double *points = mutable_points();
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < num_points_; ++j) {
            tmp[j] = points[3 * j + i];
        }
        median(i) = Median(&tmp);
    }

    for (int i = 0; i < num_points_; ++i) {
        VectorRef point(points + 3 * i, 3);
        tmp[i] = (point - median).lpNorm<1>();
    }

    const double median_absolute_deviation = Median(&tmp);

    // Scale so that the median absolute deviation of the resulting
    // reconstruction is 100

    const double scale = 100.0 / median_absolute_deviation;

    // X = scale * (X - median)
    for (int i = 0; i < num_points_; ++i) {
        VectorRef point(points + 3 * i, 3);
        point = scale * (point - median);
    }

    double *cameras = mutable_cameras();
    double angle_axis[3];
    double center[3];
    for (int i = 0; i < num_cameras_; ++i) {
        double *camera = cameras + camera_block_size() * i;
        CameraToAngelAxisAndCenter(camera, angle_axis, center);
        // center = scale * (center - median)
        VectorRef(center, 3) = scale * (VectorRef(center, 3) - median);
        AngleAxisAndCenterToCamera(angle_axis, center, camera);
    }
}

void BALProblem::Perturb(const double rotation_sigma,
                         const double translation_sigma,
                         const double point_sigma) {
    assert(point_sigma >= 0.0);
    assert(rotation_sigma >= 0.0);
    assert(translation_sigma >= 0.0);

    double *points = mutable_points();
    if (point_sigma > 0) {
        for (int i = 0; i < num_points_; ++i) {
            PerturbPoint3(point_sigma, points + 3 * i);
        }
    }

    for (int i = 0; i < num_cameras_; ++i) {
        double *camera = mutable_cameras() + camera_block_size() * i;

        double angle_axis[3];
        double center[3];
        // Perturb in the rotation of the camera in the angle-axis
        // representation
        CameraToAngelAxisAndCenter(camera, angle_axis, center);
        if (rotation_sigma > 0.0) {
            PerturbPoint3(rotation_sigma, angle_axis);
        }
        AngleAxisAndCenterToCamera(angle_axis, center, camera);

        if (translation_sigma > 0.0)
            PerturbPoint3(translation_sigma, camera + camera_block_size() - 6);
    }
}
