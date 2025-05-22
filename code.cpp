#include "code.h"

using namespace Eigen;

int Worker::i_k(Workspace& ws_) const
{
	ws_.cp = ws_.cpos_;
	ws_.cp(0) *= -1.0;

	ws_.cp += params_.N;

	ws_.ikin_10_ = ws_.cp(0) - params_.Z(0);
	ws_.ikin_11_ = ws_.cp(1) - params_.Z(1);

	ws_.ZK_ = sqrt(ws_.ikin_10_ * ws_.ikin_10_ + ws_.ikin_11_ * ws_.ikin_11_);

	ws_.ikin_1_ = 1.0 / ws_.ZK_;
	ws_.cosBeta_1_ = ws_.ikin_10_ * ws_.ikin_1_;
	ws_.sinBeta_1_ = ws_.ikin_11_ * ws_.ikin_1_;

	ws_.ikin_5_ = params_.ZCsqr - params_.c_t->r * params_.c_t->r;
	ws_.ikin_7_ = 1.0 / (2.0 * params_.ZC * ws_.ZK_);
	ws_.cosBeta_2_ = (ws_.ZK_ * ws_.ZK_ + ws_.ikin_5_) * ws_.ikin_7_;
	ws_.sinBeta_2_ = sqrt(1 - ws_.cosBeta_2_ * ws_.cosBeta_2_);

	double cosBeta = ws_.cosBeta_2_ * ws_.cosBeta_1_ - ws_.sinBeta_2_ * ws_.sinBeta_1_;
	double sinBeta = ws_.sinBeta_2_ * ws_.cosBeta_1_ + ws_.cosBeta_2_ * ws_.sinBeta_1_;

	ws_.C_(0) = params_.Z(0) + params_.ZC * cosBeta;
	ws_.C_(1) = params_.Z(1) + params_.ZC * sinBeta;

	ws_.ikin_14_ = params_.ZT * params_.ZS;
	double cosTheta = cosBeta * params_.cosTheta_s_z + sinBeta * params_.sinTheta_s_z;
	ws_.jp(0) = sqrt(params_.ZTsqr + params_.ZSsqr - 2.0 * ws_.ikin_14_ * cosTheta);

	ws_.ikin_8_ = 1.0 / params_.c_t->r;
	double cosGamma_t = (ws_.cp(0) - ws_.C_(0)) * ws_.ikin_8_;
	double sinGamma_t = (ws_.cp(1) - ws_.C_(1)) * ws_.ikin_8_;

	double cosGamma_1 = cosGamma_t * params_.c_t->cosPhi - sinGamma_t * params_.c_t->sinPhi;
	double sinGamma_1 = sinGamma_t * params_.c_t->cosPhi + cosGamma_t * params_.c_t->sinPhi;

	ws_.D_(0) = ws_.C_(0) + params_.CD * cosGamma_1;
	ws_.D_(1) = ws_.C_(1) + params_.CD * sinGamma_1;

	ws_.ikin_2_ = params_.F(0) - ws_.D_(0);
	ws_.ikin_3_ = params_.F(1) - ws_.D_(1);
	ws_.jp(1) = sqrt(ws_.ikin_2_ * ws_.ikin_2_ + ws_.ikin_3_ * ws_.ikin_3_);

	ws_.jpos_.head<2>() = ws_.jp - params_.jpos_offset;

	return 0;
}

int Worker::i_kd_using_eigen_operation(Workspace& ws_) const noexcept
{

	ws_.cvel_(0) *= -1.0;

	ws_.dZK_(0) = ws_.cosBeta_1_;
	ws_.dZK_(1) = ws_.sinBeta_1_;

	ws_.ikin_4_ = ws_.ikin_1_ * ws_.ikin_1_;
	ws_.dcosBeta_1_ =
		(ws_.ZK_ * Matrix2d::Identity().col(0) - ws_.ikin_10_ * ws_.dZK_) * ws_.ikin_4_;
	ws_.dsinBeta_1_ =
		(ws_.ZK_ * Matrix2d::Identity().col(1) - ws_.ikin_11_ * ws_.dZK_) * ws_.ikin_4_;

	ws_.ikin_9_ = 1.0 / ws_.sinBeta_2_;
	ws_.dcosBeta_2_ =
		(ws_.dZK_ *
			(ws_.ZK_ * ws_.ZK_ + params_.c_t->r * params_.c_t->r - params_.ZCsqr)) *
		ws_.ikin_7_ * ws_.ikin_1_;
	ws_.dsinBeta_2_ = -ws_.cosBeta_2_ * ws_.dcosBeta_2_ * ws_.ikin_9_;

	Eigen::Vector2d dcosBeta_ = ws_.dcosBeta_2_ * ws_.cosBeta_1_ +
		ws_.dcosBeta_1_ * ws_.cosBeta_2_ - ws_.dsinBeta_2_ * ws_.sinBeta_1_ -
		ws_.dsinBeta_1_ * ws_.sinBeta_2_;
	Eigen::Vector2d dsinBeta_ = ws_.dsinBeta_2_ * ws_.cosBeta_1_ +
		ws_.dcosBeta_1_ * ws_.sinBeta_2_ + ws_.dcosBeta_2_ * ws_.sinBeta_1_ +
		ws_.dsinBeta_1_ * ws_.cosBeta_2_;

	Eigen::Vector2d dCx = params_.ZC * dcosBeta_;
	Eigen::Vector2d dCy = params_.ZC * dsinBeta_;

	Eigen::Vector2d dcosTheta = dcosBeta_ * params_.cosTheta_s_z + dsinBeta_ *
		params_.sinTheta_s_z;

	Eigen::Vector2d dcosGamma_t = (Matrix2d::Identity().col(0) - dCx) * ws_.ikin_8_;
	Eigen::Vector2d dsinGamma_t = (Matrix2d::Identity().col(1) - dCy) * ws_.ikin_8_;

	Eigen::Vector2d dcosGamma_1 =
		dcosGamma_t * params_.c_t->cosPhi - dsinGamma_t * params_.c_t->sinPhi;
	Eigen::Vector2d dsinGamma_1 =
		dsinGamma_t * params_.c_t->cosPhi + dcosGamma_t * params_.c_t->sinPhi;

	ws_.dDx_ = dCx + params_.CD * dcosGamma_1;
	ws_.dDy_ = dCy + params_.CD * dsinGamma_1;

	ws_.jacob_.row(0) = -ws_.ikin_14_ * dcosTheta / ws_.jp(0);
	ws_.jacob_.row(1) =
		-((ws_.dDx_ * ws_.ikin_2_) + (ws_.dDy_ * ws_.ikin_3_)) / ws_.jp(1);

	ws_.jvel_.head<2>() = ws_.jacob_ * ws_.cvel_;

	return 0;
}

int Worker::i_kd_no_eigen_operation(Workspace& ws_) const noexcept
{
	ws_.cvel_(0) *= -1.0;

	ws_.dZK_(0) = ws_.cosBeta_1_;
	ws_.dZK_(1) = ws_.sinBeta_1_;

	ws_.ikin_4_ = ws_.ikin_1_ * ws_.ikin_1_;
	ws_.dcosBeta_1_(0) = (ws_.ZK_ - ws_.ikin_10_ * ws_.dZK_(0)) * ws_.ikin_4_;
	ws_.dcosBeta_1_(1) = (-ws_.ikin_10_ * ws_.dZK_(1)) * ws_.ikin_4_;

	ws_.dsinBeta_1_(0) = (-ws_.ikin_11_ * ws_.dZK_(0)) * ws_.ikin_4_;
	ws_.dsinBeta_1_(1) = (ws_.ZK_ - ws_.ikin_11_ * ws_.dZK_(1)) * ws_.ikin_4_;

	ws_.ikin_9_ = 1.0 / ws_.sinBeta_2_;
	ws_.dcosBeta_2_(0) =
		(ws_.dZK_(0) *
			(ws_.ZK_ * ws_.ZK_ + params_.c_t->r * params_.c_t->r - params_.ZCsqr)) *
		ws_.ikin_7_ * ws_.ikin_1_;
	ws_.dcosBeta_2_(1) =
		(ws_.dZK_(1) *
			(ws_.ZK_ * ws_.ZK_ + params_.c_t->r * params_.c_t->r - params_.ZCsqr)) *
		ws_.ikin_7_ * ws_.ikin_1_;

	ws_.dsinBeta_2_(0) = -ws_.cosBeta_2_ * ws_.dcosBeta_2_(0) * ws_.ikin_9_;
	ws_.dsinBeta_2_(1) = -ws_.cosBeta_2_ * ws_.dcosBeta_2_(1) * ws_.ikin_9_;

	Eigen::Vector2d dcosBeta_;
	dcosBeta_(0) = ws_.dcosBeta_2_(0) * ws_.cosBeta_1_ + ws_.dcosBeta_1_(0) * ws_.cosBeta_2_ -
		ws_.dsinBeta_2_(0) * ws_.sinBeta_1_ - ws_.dsinBeta_1_(0) * ws_.sinBeta_2_;
	dcosBeta_(1) = ws_.dcosBeta_2_(1) * ws_.cosBeta_1_ + ws_.dcosBeta_1_(1) * ws_.cosBeta_2_ -
		ws_.dsinBeta_2_(1) * ws_.sinBeta_1_ - ws_.dsinBeta_1_(1) * ws_.sinBeta_2_;

	Eigen::Vector2d dsinBeta_;
	dsinBeta_(0) = ws_.dsinBeta_2_(0) * ws_.cosBeta_1_ + ws_.dcosBeta_1_(0) * ws_.sinBeta_2_ +
		ws_.dcosBeta_2_(0) * ws_.sinBeta_1_ + ws_.dsinBeta_1_(0) * ws_.cosBeta_2_;
	dsinBeta_(1) = ws_.dsinBeta_2_(1) * ws_.cosBeta_1_ + ws_.dcosBeta_1_(1) * ws_.sinBeta_2_ +
		ws_.dcosBeta_2_(1) * ws_.sinBeta_1_ + ws_.dsinBeta_1_(1) * ws_.cosBeta_2_;

	Eigen::Vector2d dCx;
	Eigen::Vector2d dCy;
	dCx(0) = params_.ZC * dcosBeta_(0);
	dCx(1) = params_.ZC * dcosBeta_(1);
	dCy(0) = params_.ZC * dsinBeta_(0);
	dCy(1) = params_.ZC * dsinBeta_(1);

	Eigen::Vector2d dcosTheta;
	dcosTheta(0) = dcosBeta_(0) * params_.cosTheta_s_z + dsinBeta_(0) * params_.sinTheta_s_z;
	dcosTheta(1) = dcosBeta_(1) * params_.cosTheta_s_z + dsinBeta_(1) * params_.sinTheta_s_z;

	Eigen::Vector2d dcosGamma_t;
	Eigen::Vector2d dsinGamma_t;
	dcosGamma_t(0) = (1 - dCx(0)) * ws_.ikin_8_;
	dcosGamma_t(1) = (-dCx(1)) * ws_.ikin_8_;
	dsinGamma_t(0) = (-dCy(0)) * ws_.ikin_8_;
	dsinGamma_t(1) = (1 - dCy(1)) * ws_.ikin_8_;

	Eigen::Vector2d dcosGamma_1;
	dcosGamma_1(0) =
		dcosGamma_t(0) * params_.c_t->cosPhi - dsinGamma_t(0) * params_.c_t->sinPhi;
	dcosGamma_1(1) =
		dcosGamma_t(1) * params_.c_t->cosPhi - dsinGamma_t(1) * params_.c_t->sinPhi;
	Eigen::Vector2d dsinGamma_1;
	dsinGamma_1(0) =
		dsinGamma_t(0) * params_.c_t->cosPhi + dcosGamma_t(0) * params_.c_t->sinPhi;
	dsinGamma_1(1) =
		dsinGamma_t(1) * params_.c_t->cosPhi + dcosGamma_t(1) * params_.c_t->sinPhi;

	ws_.dDx_(0) = dCx(0) + params_.CD * dcosGamma_1(0);
	ws_.dDx_(1) = dCx(1) + params_.CD * dcosGamma_1(1);
	ws_.dDy_(0) = dCy(0) + params_.CD * dsinGamma_1(0);
	ws_.dDy_(1) = dCy(1) + params_.CD * dsinGamma_1(1);

	ws_.jacob_(0, 0) = -ws_.ikin_14_ * dcosTheta(0) / ws_.jp(0);
	ws_.jacob_(1, 0) = -ws_.ikin_14_ * dcosTheta(1) / ws_.jp(0);
	ws_.jacob_(1, 0) =
		-((ws_.dDx_(0) * ws_.ikin_2_) + (ws_.dDy_(0) * ws_.ikin_3_)) / ws_.jp(1);
	ws_.jacob_(1, 1) =
		-((ws_.dDx_(1) * ws_.ikin_2_) + (ws_.dDy_(1) * ws_.ikin_3_)) / ws_.jp(1);

	ws_.jvel_.head<2>() = ws_.jacob_ * ws_.cvel_;

	return 0;
}
