#ifdef __cplusplus
extern "C" {
#endif 

#include <Python/Python.h>

	/*
	 *  * Function to be called from Python
	 *   */
	static PyObject* py_myFunction(PyObject* self, PyObject* args)
	{
		    char *s = "Hello from C!";
		        return Py_BuildValue("s", s);
	}

	/*
	 *  * Bind Python function names to our C functions
	 *   */
	static PyMethodDef myModule_methods[] = {
		    {"myFunction", py_myFunction, METH_VARARGS},
		        {NULL, NULL}
	};

	/*
	 *  * Python calls this to let us initialize our module
	 *   */
	void initmyModule()
	{
		    (void) Py_InitModule("myModule", myModule_methods);
	}

#ifdef __cplusplus
}  // extern "C"
#endif 
