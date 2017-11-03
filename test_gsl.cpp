/*
 * =====================================================================================
 *
 *       Filename:  test_gsl.cpp
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  2016年11月09日 23时59分44秒
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Duan Junming (DSEC), duanjm@pku.edu.cn
 *   Organization:  
 *
 * =====================================================================================
 */

#include <gsl/gsl_integration.h>
#include <iostream>

double f(double x, void* params) {
	return *(double*)params*x*x;
}

int main() {
	gsl_function g;
	g.function = f;
	double a = 3;
	g.params = &a;
	double r, er;
	gsl_integration_workspace *w = gsl_integration_workspace_alloc(10);
	gsl_integration_qag(&g, 1, 2,	1e-6, 1e-6, 10,	1, w, &er, &r);
	gsl_integration_workspace_free(w);
	std::cout << "er: " << er << "  r: " << r << std::endl;

	return 0;
}

