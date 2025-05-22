#pragma once

#include "eigen-master/Eigen/Core"

class Workspace {
	friend class Worker;

public:
	using JointsD = Eigen::Matrix<double, 2, 1>;
	Eigen::Vector2d jp;
	Eigen::Vector2d cp;
	Eigen::Vector2d C_;
	Eigen::Vector2d dC_;
	double ZK_;
	double cosBeta_1_;
	double sinBeta_1_;
	double cosBeta_2_;
	double sinBeta_2_;
	Eigen::Vector2d D_;
	Eigen::Vector2d dZK_;
	Eigen::Vector2d dcosBeta_1_;
	Eigen::Vector2d dsinBeta_1_;
	Eigen::Vector2d dcosBeta_2_;
	Eigen::Vector2d dsinBeta_2_;
	Eigen::Vector2d dDx_;
	Eigen::Vector2d dDy_;
	double ikin_1_;
	double ikin_2_;
	double ikin_3_;
	double ikin_4_;
	double ikin_5_;
	double ikin_7_;
	double ikin_8_;
	double ikin_9_;
	double ikin_10_;
	double ikin_11_;
	double ikin_14_;
	Eigen::Matrix2d jacob_;
	JointsD jpos_;
	JointsD jvel_;
	JointsD cpos_;
	JointsD cvel_;

public:
	/** Constructor */
	Workspace()
		: jpos_(JointsD::Zero())
		, jvel_(JointsD::Zero())
		, cpos_(JointsD::Zero())
		, cvel_(JointsD::Zero())
	{
	}

	void init(JointsD& jpos, JointsD& jvel, Eigen::Vector2d& cpos,
		Eigen::Vector2d& cvel) noexcept
	{
		jpos_ = jpos;
		jvel_ = jvel;
		cpos_ = cpos;
		cvel_ = cvel;
	}
};

class Worker {

public:
	static constexpr double deg2rad = 3.14159 / 180.0;
	static constexpr double rad2deg = 180.0 / 3.14159;


private:
	struct Params {
		struct Tool {
			double r;
			double phi;
			double cosPhi;
			double sinPhi;

			void init() noexcept
			{
				cosPhi = cos(phi);
				sinPhi = sin(phi);
			}
		};
		Eigen::Vector2d jpos_offset;
		Eigen::Vector2d Z;
		Eigen::Vector2d N;
		Eigen::Vector2d S;
		Eigen::Vector2d F;
		double ZT;
		double ZC;
		double CT;
		double CD;
		double theta_z;
		double theta_s;
		Tool l_t;
		Tool u_t;
		Tool* c_t;
		double ZS;
		double ZSsqr;
		double ZTsqr;
		double ZCsqr;
		double CTsqr;
		double CDsqr;
		double cosTheta_s_z;
		double sinTheta_s_z;
		bool initialized;
		/** Constructor */
		Params() noexcept // TODO : some members are not initialized!!!
			: jpos_offset(Eigen::Vector2d::Zero())
			, Z(Eigen::Vector2d::Zero())
			, N(Eigen::Vector2d::Zero())
			, S(Eigen::Vector2d::Zero())
			, F(Eigen::Vector2d::Zero())
			, ZT(0.0)
			, ZC(0.0)
			, CT(0.0)
			, CD(0.0)
			, theta_z(0.0)
			, theta_s(0.0)
			, l_t({ 0.0, 0.0, 0.0, 0.0 })
			, u_t({ 0.0, 0.0, 0.0, 0.0 })
			, c_t(&l_t)
			, initialized(false)
		{
		}

		int init() noexcept
		{
			double tmp1 = S(0) - Z(0);
			double tmp2 = S(1) - Z(1);
			ZSsqr = tmp1 * tmp1 + tmp2 * tmp2;
			ZS = sqrt(ZSsqr);
			ZTsqr = ZT * ZT;
			ZCsqr = ZC * ZC;
			CDsqr = CD * CD;
			CTsqr = CT * CT;
			theta_s = atan2(tmp2, tmp1);
			if (theta_z == 0.0)
			{
				tmp1 = (ZTsqr + ZCsqr - CTsqr) / (2 * ZT * ZC);
				tmp2 = sqrt(1 - tmp1 * tmp1);
				theta_z = atan2(tmp2, tmp1);
			}
			cosTheta_s_z = cos(theta_s - theta_z);
			sinTheta_s_z = sin(theta_s - theta_z);
			u_t.init();
			l_t.init();
			initialized = true;

			return 0;
		}
	};

public:
	Params params_;

	Worker()
		: params_()
	{
	}


public:
	int i_kd_no_eigen_operation(Workspace& ws_) const noexcept;
	int i_kd_using_eigen_operation(Workspace& ws_) const noexcept;
	int i_k(Workspace& ws_) const;
};


class BenchmarkTEST {
public:
	Worker kp_;
	Workspace w;

	Workspace::JointsD jjpos_;
	Workspace::JointsD jjvel_;
	Eigen::Vector2d ccpos_;
	Eigen::Vector2d ccvel_;


public:
	Worker& k_;
	Workspace* ws_;


	BenchmarkTEST() noexcept
		: k_(kp_)
		, ws_(&w)
	{
	}

public:
	int setup() noexcept
	{
		if (ws_ == nullptr) return -1;
		ws_->init(jjpos_, jjvel_, ccpos_, ccvel_);

		ws_->cpos_ << 43.723, 247.101;
		ws_->cvel_ << 26.237, -76.038;

		k_.params_.jpos_offset(0) = 253.0;
		k_.params_.jpos_offset(1) = 687.0;
		k_.params_.N(0) = 1534.0;
		k_.params_.N(1) = 410.0;
		k_.params_.S(0) = 283.0;
		k_.params_.S(1) = 1100.0;
		k_.params_.F(0) = 900.0;
		k_.params_.F(1) = 1651.0;
		k_.params_.ZT = 1058.12439;
		k_.params_.ZC = 410.0;
		k_.params_.CT = 650.0;
		k_.params_.CD = 1176.542;
		k_.params_.l_t.r = 1513.311786;
		k_.params_.l_t.phi = 27.773039 * Worker::deg2rad;
		k_.params_.u_t.r = 1513.311786;
		k_.params_.u_t.phi = 27.773039 * Worker::deg2rad;
		k_.params_.c_t = &k_.params_.u_t;
		k_.params_.init();

		k_.i_k(*ws_);

		return 0;
	}

	int body_no_eigen_operation() noexcept
	{
		return k_.i_kd_no_eigen_operation(*ws_);
	}

	int body_with_eigen_operation() noexcept
	{
		return k_.i_kd_using_eigen_operation(*ws_);
	}
};