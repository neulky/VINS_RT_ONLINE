#include "initial_alignment.h"

void solveGyroscopeBias(map<double, ImageFrame> &all_image_frame, Vector3d* Bgs)
{
    Matrix3d A;
    Vector3d b;
    Vector3d delta_bg;
    A.setZero();
    b.setZero();
    map<double, ImageFrame>::iterator frame_i;
    map<double, ImageFrame>::iterator frame_j;
    for (frame_i = all_image_frame.begin(); next(frame_i) != all_image_frame.end(); frame_i++)
    {
        frame_j = next(frame_i);
        MatrixXd tmp_A(3, 3);
        tmp_A.setZero();
        VectorXd tmp_b(3);
        tmp_b.setZero();
        Eigen::Quaterniond q_ij(frame_i->second.R.transpose() * frame_j->second.R);
        tmp_A = frame_j->second.pre_integration->jacobian.template block<3, 3>(O_R, O_BG);
        tmp_b = 2 * (frame_j->second.pre_integration->delta_q.inverse() * q_ij).vec();
        A += tmp_A.transpose() * tmp_A;
        b += tmp_A.transpose() * tmp_b;

    }
    delta_bg = A.ldlt().solve(b);
    ROS_WARN_STREAM("gyroscope bias initial calibration " << delta_bg.transpose());

    for (int i = 0; i <= WINDOW_SIZE; i++)
        Bgs[i] += delta_bg;

    for (frame_i = all_image_frame.begin(); next(frame_i) != all_image_frame.end( ); frame_i++)
    {
        frame_j = next(frame_i);
        frame_j->second.pre_integration->repropagate(Vector3d::Zero(), Bgs[0]);
    }
}


MatrixXd TangentBasis(Vector3d &g0)
{
    Vector3d b, c;
    Vector3d a = g0.normalized();
    Vector3d tmp(0, 0, 1);
    if(a == tmp)
        tmp << 1, 0, 0;
    b = (tmp - a * (a.transpose() * tmp)).normalized();
    c = a.cross(b);
    MatrixXd bc(3, 2);
    bc.block<3, 1>(0, 0) = b;
    bc.block<3, 1>(0, 1) = c;
    return bc;
}

void RefineGravity(map<double, ImageFrame> &all_image_frame, Vector3d &g, VectorXd &x)
{
    Vector3d g0 = g.normalized() * G.norm();
    Vector3d lx, ly;
    //VectorXd x;
    int all_frame_count = all_image_frame.size();
    int n_state = all_frame_count * 3 + 2 + 1;

    MatrixXd A{n_state, n_state};
    A.setZero();
    VectorXd b{n_state};
    b.setZero();

    map<double, ImageFrame>::iterator frame_i;
    map<double, ImageFrame>::iterator frame_j;
    for(int k = 0; k < 4; k++)
    {
        MatrixXd lxly(3, 2);
        lxly = TangentBasis(g0);
        int i = 0;
        for (frame_i = all_image_frame.begin(); next(frame_i) != all_image_frame.end(); frame_i++, i++)
        {
            frame_j = next(frame_i);

            MatrixXd tmp_A(6, 9);
            tmp_A.setZero();
            VectorXd tmp_b(6);
            tmp_b.setZero();

            double dt = frame_j->second.pre_integration->sum_dt;


            tmp_A.block<3, 3>(0, 0) = -dt * Matrix3d::Identity();
            tmp_A.block<3, 2>(0, 6) = frame_i->second.R.transpose() * dt * dt / 2 * Matrix3d::Identity() * lxly;
            tmp_A.block<3, 1>(0, 8) = frame_i->second.R.transpose() * (frame_j->second.T - frame_i->second.T) / 100.0;     
            tmp_b.block<3, 1>(0, 0) = frame_j->second.pre_integration->delta_p + frame_i->second.R.transpose() * frame_j->second.R * TIC[0] - TIC[0] - frame_i->second.R.transpose() * dt * dt / 2 * g0;

            tmp_A.block<3, 3>(3, 0) = -Matrix3d::Identity();
            tmp_A.block<3, 3>(3, 3) = frame_i->second.R.transpose() * frame_j->second.R;
            tmp_A.block<3, 2>(3, 6) = frame_i->second.R.transpose() * dt * Matrix3d::Identity() * lxly;
            tmp_b.block<3, 1>(3, 0) = frame_j->second.pre_integration->delta_v - frame_i->second.R.transpose() * dt * Matrix3d::Identity() * g0;


            Matrix<double, 6, 6> cov_inv = Matrix<double, 6, 6>::Zero();
            //cov.block<6, 6>(0, 0) = IMU_cov[i + 1];
            //MatrixXd cov_inv = cov.inverse();
            cov_inv.setIdentity();

            MatrixXd r_A = tmp_A.transpose() * cov_inv * tmp_A;
            VectorXd r_b = tmp_A.transpose() * cov_inv * tmp_b;

            A.block<6, 6>(i * 3, i * 3) += r_A.topLeftCorner<6, 6>();
            b.segment<6>(i * 3) += r_b.head<6>();

            A.bottomRightCorner<3, 3>() += r_A.bottomRightCorner<3, 3>();
            b.tail<3>() += r_b.tail<3>();

            A.block<6, 3>(i * 3, n_state - 3) += r_A.topRightCorner<6, 3>();
            A.block<3, 6>(n_state - 3, i * 3) += r_A.bottomLeftCorner<3, 6>();
        }
            A = A * 1000.0;
            b = b * 1000.0;
            x = A.ldlt().solve(b);
            VectorXd dg = x.segment<2>(n_state - 3);
            g0 = (g0 + lxly * dg).normalized() * G.norm();
            //double s = x(n_state - 1);
    }   
    g = g0;
}

/* bool LinearAlignment(map<double, ImageFrame> &all_image_frame, Vector3d &g, VectorXd &x)
{
    int all_frame_count = all_image_frame.size();
    int n_state = all_frame_count * 3 + 3 + 1;

    MatrixXd A{n_state, n_state};
    A.setZero();
    VectorXd b{n_state};
    b.setZero();

    map<double, ImageFrame>::iterator frame_i;
    map<double, ImageFrame>::iterator frame_j;
    int i = 0;
    for (frame_i = all_image_frame.begin(); next(frame_i) != all_image_frame.end(); frame_i++, i++)
    {
        frame_j = next(frame_i);

        MatrixXd tmp_A(6, 10);
        tmp_A.setZero();
        VectorXd tmp_b(6);
        tmp_b.setZero();

        double dt = frame_j->second.pre_integration->sum_dt;

        tmp_A.block<3, 3>(0, 0) = -dt * Matrix3d::Identity();
        tmp_A.block<3, 3>(0, 6) = frame_i->second.R.transpose() * dt * dt / 2 * Matrix3d::Identity();
        tmp_A.block<3, 1>(0, 9) = frame_i->second.R.transpose() * (frame_j->second.T - frame_i->second.T) / 100.0;     
        tmp_b.block<3, 1>(0, 0) = frame_j->second.pre_integration->delta_p + frame_i->second.R.transpose() * frame_j->second.R * TIC[0] - TIC[0];
        //cout << "delta_p   " << frame_j->second.pre_integration->delta_p.transpose() << endl;
        tmp_A.block<3, 3>(3, 0) = -Matrix3d::Identity();
        tmp_A.block<3, 3>(3, 3) = frame_i->second.R.transpose() * frame_j->second.R;
        tmp_A.block<3, 3>(3, 6) = frame_i->second.R.transpose() * dt * Matrix3d::Identity();
        tmp_b.block<3, 1>(3, 0) = frame_j->second.pre_integration->delta_v;
        //cout << "delta_v   " << frame_j->second.pre_integration->delta_v.transpose() << endl;

        Matrix<double, 6, 6> cov_inv = Matrix<double, 6, 6>::Zero();
        //cov.block<6, 6>(0, 0) = IMU_cov[i + 1];
        //MatrixXd cov_inv = cov.inverse();
        cov_inv.setIdentity();

        MatrixXd r_A = tmp_A.transpose() * cov_inv * tmp_A;
        VectorXd r_b = tmp_A.transpose() * cov_inv * tmp_b;

        A.block<6, 6>(i * 3, i * 3) += r_A.topLeftCorner<6, 6>();
        b.segment<6>(i * 3) += r_b.head<6>();

        A.bottomRightCorner<4, 4>() += r_A.bottomRightCorner<4, 4>();
        b.tail<4>() += r_b.tail<4>();

        A.block<6, 4>(i * 3, n_state - 4) += r_A.topRightCorner<6, 4>();
        A.block<4, 6>(n_state - 4, i * 3) += r_A.bottomLeftCorner<4, 6>();
    }
    A = A * 1000.0;
    b = b * 1000.0;
    x = A.ldlt().solve(b);
    double s = x(n_state - 1) / 100.0;
    ROS_DEBUG("estimated scale: %f", s);
    g = x.segment<3>(n_state - 4);
    ROS_DEBUG_STREAM(" result g     " << g.norm() << " " << g.transpose());
    if(fabs(g.norm() - G.norm()) > 1.0 || s < 0)
    {
        return false;
    }

    RefineGravity(all_image_frame, g, x);
    s = (x.tail<1>())(0) / 100.0;
    (x.tail<1>())(0) = s;
    ROS_DEBUG_STREAM(" refine     " << g.norm() << " " << g.transpose());
    if(s < 0.0 )
        return false;   
    else
        return true;
} */

bool LinearAlignment(map<double, ImageFrame> &all_image_frame, Vector3d &g, double &s, Vector3d &ba)
{
    cv::Mat Rbc = Converter::toCvMat(TIC[0]);
    cv::Mat Rcb = Rbc.t();
    
    //Step 1.
    //Approx Scale and Gravity vector in 'world' frame (first KF's camera frame)
    int frame_count = all_image_frame.size();
    //Solve A*x=B for x=[s,gw,Pcb]
    cv::Mat A = cv::Mat::zeros(3*(frame_count-2), 7, CV_64F);
    cv::Mat B = cv::Mat::zeros(3*(frame_count-2), 1, CV_64F);
    cv::Mat I3 = cv::Mat::eye(3, 3, CV_64F);

    map<double, ImageFrame>::iterator frame_i;
    map<double, ImageFrame>::iterator frame_j;
    map<double, ImageFrame>::iterator frame_k;
    int i = 0;
    for(frame_i = all_image_frame.begin(); next(next(frame_i)) != all_image_frame.end();frame_i++,i++)
    {
        frame_j = next(frame_i);
        frame_k = next(frame_j);
        //Delta time between frames
        double dt12 = frame_j->second.pre_integration->sum_dt;
        double dt23 = frame_k->second.pre_integration->sum_dt;
        //Pre-integrated measurements
        cv::Mat dp12 = Converter::toCvMat(frame_j->second.pre_integration->delat_p);
        cv::Mat dv12 = Converter::toCvMat(frame_j->second.pre_integration->delat_v);
        cv::Mat dp23 = Converter::toCvMat(frame_k->second.pre_integration->delat_p);
        //Position of camera
        cv::Mat pc1 = frame_i->second.T;
        cv::Mat pc2 = frame_j->second.T;
        cv::Mat pc3 = frame_k->second.T;
        //Rotation of camera
        cv::Mat Rc1 = Converter::toCvMat(frame_i->second.R * RIC[0]);
        cv::Mat Rc2 = Converter::toCvMat(frame_j->second.R * RIC[0]);
        cv::Mat Rc3 = Converter::toCvMat(frame_k->second.R * RIC[0]);

        //Stack to A/B matrix
        //lambda*s + beta*g + erfa*Pcb = gamma
        cv::Mat lambda = (pc2 - pc1) * dt23 + (pc2 - pc3) * dt12;
        cv::Mat beta = 0.5 * I3 * (dt12 * dt12 * dt23) + dt12 * dt23 * dt23);
        cv::Mat erfa = (Rc2 - Rc3) * dt12 - (Rc1 - Rc2) * dt23;
        cv::Mat gamma = Rc1 * Rcb * dp12 * dt23 - Rc2 * Rcb * dp23 * dt12 - Rc1 * Rcb * dv12 * dt12 * dt23;
        lambda.copyTo(A.rowRang(3*i+0,3*i+3).col(0));
        beta.copyTo(A.rowRange(3*i+0,3*i+3).colRange(1,4));
        erfa.copyTo(A.rowRange(3*i+0,3*i+3).colRange(4,7));
        gamma.copyTo(B.rowRange(3*i+0,3*i+3));
    }
    //Use svd to compute A*x=B, x=[s,gw,Pcb] 7x1 vector
    //A = u*w*vt, u*w*vt*x=B
    //Then x = vt'*winv*u'*B
    cv::Mat w, u, vt;
    //Ais changed in SVDecomp() with cv::SVD::MODIFY_A for speed
    cv::SVDecomp(A, w, u, vt, cv::SVD::MODIFY_A);
    //Compute winv
    cv::Mat winv = cv::Mat::eye(4, 4, CV_64F);
    for(int i = 0;i < 7;i++)
    {
        if(fabs(w.at<double>(i)) < 1e-10)
        {
            w.at<double>(i) += 1e-10;
        }
        winv.at<double>(i, i) = 1./w.at<double>(i);
    }
    //Then x = vt'*winv*u'*B
    cv::Mat x = vt.t() * winv * u.t() * B;

    //x = [s, gw, Pcb] 7x1 vector
    double sstar = x.at<double>(0);
    cv::Mat gwstar = x.rowRange(1, 4);
    cv::Mat Pcb1 = x.rowRange(4, 7);

    //Step 3.
    //Use gravity magnitude 9.8 as constraint
    //gI = [0,0,1], the normalized gravity vector in an inertial frame, NED type with no orientation.
    cv::Mat gI = cv::Mat::zeros(3, 1, CV_64F);
    gI.at<double>(2) = 1;
    //Normalized approx. gravity vector in world frame
    cv::Mat gwn = gwstar/cv::norm(gwstar);
    //vhat = (gI x gw)/|gI x gw|
    cv::Mat gIxgwn = gI.cross(gwn);
    double normgIxgwn = cv::norm(gIxgwn);
    cv::Mat vhat = gIxgwn/normgIxgwn;
    double theta = std::atan2(normgIxgwn, gI.dot(gwn));

    //RwI
    Eigen::Vector3d vhateig = Converter::toVector3d(vhat);
    Eigen::Matrix3d RWIeig = Sophus::SO3::exp(vhateig*theta).matrix();
    cv::Mat Rwi = Converter::toCvMat(RWIeig);
    cv::Mat GI = gI * 9.8012;
    //Solve C*x=D for x = [s, dthetaxy, ba, Pcb] 9x1 vector
    cv::Mat C = cv::Mat::zeros(3*(frame_count-2), 9, CV_64F);
    cv::Mat D = cv::Mat::zeros(3*(frame_count-2), 1, CV_64F);

    i = 0;
    for(frame_i = all_image_frame.begin();next(next(frame_i)) != all_image_frame.end();frame_i++,i++))
    {
        //Delta time between frames
        double dt12 = frame_j->second.pre_integration->sum_dt;
        double dt23 = frame_k->second.pre_integration->sum_dt;
        //Pre-integrated measurements
        Eigen::Matrix3d dp_dba12 = frame_j->second.pre_integration->jacobian.block<3, 3>(O_P, O_BA);
        Eigen::Matrix3d dp_dba23 = frame_k->second.pre_integration->jacobian.block<3, 3>(O_P, O_BA);
        Eigen::Matrix3d dv_dba12 = frame_j->second.pre_integration->jacobian.block<3, 3>(O_V, O_BA);
        
        cv::Mat dp12 = Converter::toCvMat(frame_j->second.pre_integration->delta_p);
        cv::Mat dv12 = Converter::toCvMat(frame_j->second.pre_integration->delta_v);
        cv::Mat dp23 = Converter::toCvMat(frame_k->second.pre_integration->delta_p);
        cv::Mat Jpba12 = Converter::toCvMat(dp_dba12);
        cv::Mat Jvba12 = Converter::toCvMat(dv_dba12);
        cv::Mat Jpba23 = Converter::toCvMat(dp_dba23);
        //Position of camera
        cv::Mat pc1 = frame_i->second.T;
        cv::Mat pc2 = frame_j->second.T;
        cv::Mat pc3 = frame_k->second.T;
        //Rotation of camera
        cv::Mat Rc1 = Converter::toCvMat(frame_i->second.R * RIC[0]);
        cv::Mat Rc2 = Converter::toCvMat(frame_j->second.R * RIC[0]);
        cv::Mat Rc3 = Converter::toCvMat(frame_k->second.R * RIC[0]);
        //Stack to C/D matrix
        //lamda * s + phi * dthetaxy + zeta * ba + erfa * Pcb = psi
        cv::Mat lambda = (pc2 - pc1) * dt23 - (pc3 - pc2) * dt12;
        cv::Mat phi = -0.5 * Rwi * SkewSymmetricMatrix(GI) * (dt12 * dt23 * dt23 + dt12 * dt12 + dt23);
        cv::Mat zeta = Rc2 * Jpba23 * dt12 + Rc1 * Rcb * Jvba23 * dt12 * dt23 - Rc1 * Jpba12 * dt23;
        cv::Mat erfa = (Rc2 - Rc3) * dt12 - (Rc1 - Rc2) * dt23;
        cv::Mat psi = Rc1 * Rcb * dp12 * dt23 - Rc2 * Rcb * dp23 * dt12 - Rc1 * Rcb * dv12 * dt12 * dt23
                     - 0.5 * Rwi * GI * (dt12 * dt23 * dt23 + dt12 * dt12 * dt23);

        lambda.copyTo(C.rowRange(3*i+0,3*i+3).col(0));
        phi.colRange(0,2).copyTo(C.rowRange(3*i+0,3*i+3).colRange(1,3));
        zeta.copyTo(C.rowRange(3*i+0,3*i+3).colRange(3,6));
        erfa.copyTo(C.rowRange(3*i+0,3*i+3).colRange(6,9));
        psi.copyTo(D.rowRange(3*i+0,3*i+3));
    }

    //Use svd to compute C*x=D, x=[s, dthetaxy, ba, Pcb] 9x1 vector
    //C = u*w*vt, u*w*vt*x=D
    //Then x = vt'*Winv*u'*D
    cv::Mat w2,u2,vt2;
    //C is changed in SVDecomp() with cv::SVD::MODIFY_A for speed
    cv::SVDecomp(C, w2, u2, vt2, cv::SVD::MODIFY_A);

    //Compute winv
    cv::Mat w2inv = cv::Mat::eye(9, 9, CV_64F);
    for(int i = 0;i < 9;i++)
    {
        if(fabs(w2.at<double>(i) < 1e-10))
        {
            w2.at<double>(i) += 1e-10;
        }
        w2inv.at<double>(i,i) = 1./w2.at<double>(i);
    }

    //Then y = vt'*winv*u'*D
    cv::Mat y = vt2.t()*w2inv*u2.t()*D;

    s = y.at<double>(0);                                                //s

    //dtheta = [dx,du,0]                                                //g     
    cv::Mat dthetaxy = y.rowRange(1, 3);
    cv::Mat dtheda = cv::Mat::zeros(3, 1, CV_64F);
    dthetaxy.copyTo(dtheda.rowRange(0, 2));
    Eigen::Vector3d dthetaeig = Converter::toVector3d(dtheta);
    //Rwi_ = Rwi*exp(dtheta)
    Eigen::Matrix3d Rwieig_ = RWIeig*Sophus::SO3::exp(dthetaeig).matrix();
    cv::Mat Rwi_ = Converter::toCvMat(Rwieig_);
    cv::Mat gw_cv = Rwi_*GI;
    g = Converter::toVector3d(gw_cv);
    
    cv::Mat dbiasa_ = y.rowRange(3, 6);                                 //ba
    Vector3d dba = Converter::toVector3d(dbiasa_);

    cv::Mat Pcb_ = y.rowRange(6, 9);                                    //Pbc
    cv::Mat Pbc_ = -Rbc * Pcb_;  
    TIC[0] = Converter::toVector3d(Pbc_);

}
bool VisualIMUAlignment(map<double, ImageFrame> &all_image_frame, Vector3d* Bgs, Vector3d &g, double &s, Vector3d &ba)
{
    solveGyroscopeBias(all_image_frame, Bgs);

    if(LinearAlignment(all_image_frame, g, s, ba))
        return true;
    else 
        return false;
}
