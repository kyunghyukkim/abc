#include <iostream>

extern "C" {
	char* f(char*);
	}

void g(char* p)
{
	
}

char* f(char* p)
{
	g(p);
	return p;
}

