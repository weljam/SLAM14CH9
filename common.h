#pragma once

/// 从文件读入BAL dataset
class BALProblem {
public:
    /// load bal data from text file
    //explicit限制不能发生隐式转换
    explicit BALProblem(const std::string &filename, bool use_quaternions = false);

    ~BALProblem() {
        delete[] point_index_;
        delete[] camera_index_;
        delete[] observations_;
        delete[] parameters_;
    }

    /// save results to text file
    void WriteToFile(const std::string &filename) const;

    /// save results to ply pointcloud
    void WriteToPLYFile(const std::string &filename) const;

    void Normalize();

    void Perturb(const double rotation_sigma,
                 const double translation_sigma,
                 const double point_sigma);

    int camera_block_size() const { return use_quaternions_ ? 10 : 9; } //返回camera维数

    int point_block_size() const { return 3; }      //返回点的维数

    int num_cameras() const { return num_cameras_; }    //返回相机数

    int num_points() const { return num_points_; }      //返回路标点个数

    int num_observations() const { return num_observations_; }      //返回观测点个数

    int num_parameters() const { return num_parameters_; }      //返回位姿和路标点的个数

    const int *point_index() const { return point_index_; }     //返回像素点对应路标点的索引

    const int *camera_index() const { return camera_index_; }   //返回像素点对应位姿的索引

    const double *observations() const { return observations_; }    //返回观测点（像素）的地址

    const double *parameters() const { return parameters_; }    //参数{pose和point}的首地址

    const double *cameras() const { return parameters_; }       //位姿首地址

    const double *points() const { return parameters_ + camera_block_size() * num_cameras_; } //路标点的首地址   

    /// camera参数的起始地址
    double *mutable_cameras() { return parameters_; }       //相机参数首地址

    double *mutable_points() { return parameters_ + camera_block_size() * num_cameras_; }  //路标点的首地址

    double *mutable_camera_for_observation(int i) {
        return mutable_cameras() + camera_index_[i] * camera_block_size();      //第i个观测值对应的相机参数地址
    }

    double *mutable_point_for_observation(int i) {
        return mutable_points() + point_index_[i] * point_block_size();         //第i个观测值对应的路标点地址
    }

    const double *camera_for_observation(int i) const {
        return cameras() + camera_index_[i] * camera_block_size();              // 第i个观测值对应的相机参数地址
    }

    const double *point_for_observation(int i) const {
        return points() + point_index_[i] * point_block_size();                  //第i个观测值对应的路标点地址
    }

private:
    //
    void CameraToAngelAxisAndCenter(const double *camera,
                                    double *angle_axis,
                                    double *center) const;

    void AngleAxisAndCenterToCamera(const double *angle_axis,
                                    const double *center,
                                    double *camera) const;

    int num_cameras_;       //相机个数
    int num_points_;        //路标点个数
    int num_observations_;  //观测点个数
    int num_parameters_;
    bool use_quaternions_;

    int *point_index_;      // 每个observation对应的point index
    int *camera_index_;     // 每个observation对应的camera index
    double *observations_;
    double *parameters_;
};
