#include <stdio.h>
#include<iostream>
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
   test = "Hello char string";
   return test;
}