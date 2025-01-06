#include "udf.h"

DEFINE_MASS_TRANSFER(liq_gas_source, cell, thread, from_index, from_species_index, to_index, to_species_index)
{
	/*定义气液相的thread*/
	Thread* gas, * liq;
	gas = THREAD_SUB_THREAD(thread, from_index);
	liq = THREAD_SUB_THREAD(thread, to_index);

	/*定义变量*/
	real k_L_O2;
	real k_L_O2_LS;
	real k_L_O2_higbie;
	real MW_O2 = 31.99;
	real MW_N2 = 28.01;
	real c_s_O2;
	real y_O2;
	real c_O2;
	real m_lg_O2;
	real a;

	/*局部条件*/
	real eG = C_VOF(cell, gas);/*获取网格中气含率，也就是气体的体积分数*/
	real P = C_P(cell, thread) + 101325; /* Pa, 反应器内的压力*/
	real dissipation = C_D(cell, liq); /* 湍流耗散率, m2/s3 */
	real w_O2 = C_YI(cell, gas, 0); /* 气相中氧气的质量分数 */
	real w_N2 = C_YI(cell, gas, 1); /* 气相中氮气的质量分数 */


	/* 液相特性 */
	double density_L = C_R(cell, liq);  /* kg/m3, 液相密度 */
	real mu_l = C_MU_L(cell, liq); /* Pa.s */
	real nu_l = mu_l / density_L; /* m2/s */

	/* 气相特性 */
	real d_b = C_PHASE_DIAMETER(cell, gas);  /* 气泡直径, m */

	/* 氧气在液相中的扩散系数 */
	double D_L_O2 = 2.42 * pow(10, -9); /* m^2/s, 氧气在液相中的扩散系数 */

	/* 亨利常数*/
	double H_O2 = 4.16 * pow(10, -7); /* kg/m^3/Pa 氧气在水中的亨利常数*/

	/* 气泡的比表面积 */
	a = 6 * eG / d_b;    /* 1/m */

	/*计算传质kl，小涡计算公式*/
	k_L_O2_LS = 0.4 * pow(D_L_O2, 0.5) * pow(dissipation / nu_l, 0.25); /* m/s */

	/* 计算传质kl，用这个滑移速度模型*/
	real u_slip = C_U(cell, gas) - C_U(cell, liq); /* m/s */
	real v_slip = C_V(cell, gas) - C_V(cell, liq); /* m/s */
	real w_slip = C_W(cell, gas) - C_W(cell, liq); /* m/s */
	real vel_slip = pow(pow(u_slip, 2) + pow(v_slip, 2) + pow(w_slip, 2), 0.5); /* m/s */
	k_L_O2_higbie = pow(4 * D_L_O2 * vel_slip / d_b / 3.14, 0.5); /* m/s */

	/* 选最大值 */
	if (k_L_O2_LS > k_L_O2_higbie)
	{
		k_L_O2 = k_L_O2_LS;
	}
	else
	{
		k_L_O2 = k_L_O2_higbie;
	}
	/* Limiting k_L if k_L > 0.01] */
	if (k_L_O2 > 0.01)
	{
		k_L_O2 = 0.01;
	}
	/* 计算浓度 */
	y_O2 = (w_O2 / MW_O2) / ((w_O2 / MW_O2) + (w_N2 / MW_N2));  /* -, O2 gas mmole fraction*/
	c_s_O2 = H_O2 * P * y_O2;               /* kg/m^3 */
	c_O2 = density_L * C_YI(cell, liq, 0);    /* kg/m^3 */
	/* 计算传质速率 */
	m_lg_O2 = k_L_O2 * a * (c_s_O2 - c_O2); /* kg/m^3/s */

	return (m_lg_O2);
}
