#include <stdio.h>
#include<iostream>
using namespace std;

extern "C" {
	void TestVoid(void);
}

void TestVoid(void)
{
	
   cout << "Hello World!";
   
}