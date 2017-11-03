/*
 * =====================================================================================
 *
 *       Filename:  system_test.cpp
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  2016年06月02日 12时49分17秒
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Duan Junming (DSEC), duanjm@pku.edu.cn
 *   Organization:  
 *
 * =====================================================================================
 */

#include <iostream>
#include <cstdlib>

int main() {
	std::cout << "1" << std::endl;
	system("mv 1.dat 2.dat");
	std::cout << "2" << std::endl;
	return 0;
}

