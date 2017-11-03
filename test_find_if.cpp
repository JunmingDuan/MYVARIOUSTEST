/*
 * =====================================================================================
 *
 *       Filename:  test.cpp
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  2016年10月24日 00时09分45秒
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Duan Junming (DSEC), duanjm@pku.edu.cn
 *   Organization:  
 *
 * =====================================================================================
 */

#include <algorithm>
#include <iostream>
#include <functional>
#include <vector>

class A  {  
	public:  
		A(const std::vector<int>& str,int id)  
		{  
			this->str=str;  
			this->id=id;  
		}  
		std::vector<int> str;  
		int id;  
};

struct compare: std::binary_function<A, std::vector<int>,bool>  
{  
	bool operator()( A &value, std::vector<int> str) const  
	{  
		if (value.str == str)  
			return true;  
		else  
			return false;  
	}  
};  

int main() {
	std::vector<A> a;  
	std::vector<int> tmp;

	tmp.push_back(1);
	tmp.push_back(2);
	A b(tmp,4);  
	tmp.push_back(3);
	tmp.push_back(4);
	A c(tmp,6);  
	tmp.push_back(5);
	tmp.push_back(6);
	A d(tmp,7);  

	a.push_back(b);  
	a.push_back(c);  
	a.push_back(d);  

	std::vector<A>::iterator t=find_if(a.begin(),a.end(),bind2nd(compare(),tmp));
	std::cout << a[t-a.begin()].id << std::endl;

	return 0;
}

