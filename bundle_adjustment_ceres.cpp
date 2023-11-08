#include <iostream>
#include <ceres/ceres.h>
#include "common.h"
#include "SnavelyReprojectionError.h"

using namespace std;

void SolveBA(BALProblem &bal_problem);

int main(int argc, char **argv) {
    if (argc != 2) {
        cout << "usage: bundle_adjustment_ceres bal_data.txt" << endl;
        return 1;
    }
    //初始化得到相机参数，路标点，观测点
    BALProblem bal_problem(argv[1]);
    //归一化
    bal_problem.Normalize();
    //添加噪声
    bal_problem.Perturb(0.1, 0.5, 0.5);
    //写道meshlab的。ply文件
    bal_problem.WriteToPLYFile("initial.ply");
    SolveBA(bal_problem);
    bal_problem.WriteToPLYFile("final.ply");

    return 0;
}

void SolveBA(BALProblem &bal_problem) {
    //路标点的内存大小
    const int point_block_size = bal_problem.point_block_size();
    //相机参数的大小
    const int camera_block_size = bal_problem.camera_block_size();
    //路标点的头地址
    double *points = bal_problem.mutable_points();
    //相机参数的头地址
    double *cameras = bal_problem.mutable_cameras();

    // Observations is 2 * num_observations long array observations
    // [u_1, u_2, ... u_n], where each u_i is two dimensional, the x
    // and y position of the observation.
    //观测值
    const double *observations = bal_problem.observations();
    ceres::Problem problem;

    for (int i = 0; i < bal_problem.num_observations(); ++i) {
        ceres::CostFunction *cost_function;

        // Each Residual block takes a point and a camera as input
        // and outputs a 2 dimensional Residual
        //相当于创建代价函数，并初始化变量
        cost_function = SnavelyReprojectionError::Create(observations[2 * i + 0], observations[2 * i + 1]);

        //添加huber核函数
        ceres::LossFunction *loss_function = new ceres::HuberLoss(1.0);

        // Each observation corresponds to a pair of a camera and a point
        // which are identified by camera_index()[i] and point_index()[i]
        // respectively.
        //下一组数据
        double *camera = cameras + camera_block_size * bal_problem.camera_index()[i];
        double *point = points + point_block_size * bal_problem.point_index()[i];
        //添加误差
        problem.AddResidualBlock(cost_function, loss_function, camera, point);
    }

    // show some information here ...
    std::cout << "bal problem file loaded..." << std::endl;
    std::cout << "bal problem have " << bal_problem.num_cameras() << " cameras and "
              << bal_problem.num_points() << " points. " << std::endl;
    std::cout << "Forming " << bal_problem.num_observations() << " observations. " << std::endl;

    std::cout << "Solving ceres BA ... " << endl;
    //配置求解器
    ceres::Solver::Options options;
    //使用SPARSE——SCHUR对稀疏性加速
    options.linear_solver_type = ceres::LinearSolverType::SPARSE_SCHUR;
    //输出到cout
    options.minimizer_progress_to_stdout = true;
    //优化信息
    ceres::Solver::Summary summary;
    //开始优化
    ceres::Solve(options, &problem, &summary);
    std::cout << summary.FullReport() << "\n";
}