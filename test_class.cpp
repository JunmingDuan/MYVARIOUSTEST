/*
 * =====================================================================================
 *
 *       Filename:  test_class.cpp
 *
 *    Description:  测试在一个类中调用另一个类的函数修改其值
 *
 *        Version:  1.0
 *        Created:  2017年01月19日 21时04分36秒
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Duan Junming (DSEC), duanjm@pku.edu.cn
 *   Organization:  
 *
 * =====================================================================================
 */

#include <iostream>

class A {
	public:
		double x;

		void set_x(double y) {
			x = y;
		}
};

class B {
	public:
		A a;

		B(double init) { a.x = init; }
		void change(double aa) { a.set_x(aa); }
};

int main() {
	B b(1);
	std::cout << b.a.x << std::endl;
	b.change(2);
	std::cout << b.a.x << std::endl;

	return 0;
}

