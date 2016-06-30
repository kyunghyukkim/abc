#include <iostream>

extern "C" {
	float f(float);
	}

float f(float p)
{
	std::cout << "Hello" << std::endl;
	std::cout << "the value of p = " << p << std::endl;
	return p;
}

