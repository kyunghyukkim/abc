#include <stdio.h>
#include<iostream>
// setprecision example
#include <iostream>     // std::cout, std::fixed
#include <iomanip>      // std::setprecision

using namespace std;

extern "C" {
	
   __declspec(dllexport) void TestVoid(void);
   __declspec(dllexport) char * TestStr(char*);
   
}


void TestVoid(void)
{

   cout << "HelloWorld\n";

}

char * TestStr(char *str)
{
	
   // cout << "str\n";
   // cout << str;
   // cout << "\n";
   // cout << "*str\n";
   // cout << *str;
   // cout << "\n";
   // cout << "Hello World!";
   
   char * test;
   test = "Hello char string\n";
   
   cout << test;
   
   float Intest = 3.141579;
   
   std::cout << std::setprecision(5) << Intest << '\n';
   
   cout << Intest;
   cout << "\n";
   cout << Intest << endl;
   cout << "\n";
   
   
   return test;
}