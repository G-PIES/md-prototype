/*
 * Program: 	molecular_dynamics_OpenCL.cpp
 * Summary: 	This program is a proof of concept for both a Molecular Dynamic model
 *				as well as a test of various concepts and techniques used with OpenCL/OpenGL
 * Programmer:	Sean B. Higgins
 * Start Date:	November 13, 2023
 * Update Log:	
 * 				November 25, 2023:
 * 				Began implementation of parallel MD model. Implementing space for atoms.
 */


// 1. Program headers

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <omp.h>

#include "cl.h"
#include "cl_platform.h"

// the matrix-width and the number of work-items per work-group:
// note: the matrices are actually NUM_ATOMS * NUM_ATOMS * NUM_ATOMS and the work group sizes are LOCALSIZExLOCALSIZE:

#ifndef NUM_ELEMS
#define NUM_ELEMS		32
#endif

#ifndef NUM_ATOMS
#define NUM_ATOMS		NUM_ELEMS * NUM_ELEMS * NUM_ELEMS
#endif

#ifndef NUM_CELLS
#define	NUM_CELLS		8000
#endif

#ifndef NUM_THREADS	
#define	NUM_THREADS		128
#endif

#ifndef NUM_MODELS
#define NUM_MODELS		1
#endif

// OpenCL objects:
cl_platform_id		Platform;
cl_device_id		Device;
cl_kernel			Kernel;
cl_program			Program;
cl_context			Context;
cl_command_queue	CmdQueue;


// Do we want to print in csv file format?

//#define CSV


float			hLatice[NUM_ATOMS][NUM_ATOMS];
float			hC[NUM_ATOMS][NUM_ATOMS];
int				atomsPerCell[NUM_CELLS];
int				atomsPerLayer[NUM_ATOMS];

// The cell-layer occupancy matrix 'H' can be created.

// The atomic neighbor matrix 'W'

const char *	CL_FILE_NAME = { "molecular_dynamics.cl" };

// function prototypes:
void			SelectOpenclDevice();
char *			Vendor( cl_uint );
char *			Type( cl_device_type );
void			Wait( cl_command_queue );
int				maxVal( int*, int);


int main( int argc, char *argv[ ] )
{
#ifndef _OPENMP
        fprintf( stderr, "OpenMP is not enabled!\n" );
        return 1;
#endif

	// see if we can even open the OpenCL kernel programs
	// (no point going on if we can't):

	FILE *fp;
#ifdef WIN32
	errno_t err = fopen_s( &fp, CL_FILE_NAME, "r" );
	if ( err != 0)
#else
	fp = fopen( CL_FILE_NAME, "r" );
	if ( fp == NULL)
#endif
	{
		fprintf( stderr, "Cannot open OpenCL source file '%s'\n", CL_FILE_NAME );
		return 1;
	}

	cl_int status;		// returned status from opencl calls -- test against CL_SUCCESS


	// Get the platform id and the device id:

	SelectOpenclDevice();		// sets the global variables Platform and Device


	// 2. Allocate the host memory buffers:
	// already done -- we did it as global variables instead of on the heap so could allocate them as a 2D array

	// Initialize the input matrix:

	for( int i = 0; i < NUM_ATOMS; i++ )
	{
		for( int j = 0; j < NUM_ATOMS; j++ )
		{
			hLatice[i][j] = 1.0;
		}
	}

	

	// Initialize the count of atoms in each cell from 0 to NUM_CELLS
	// Likely can be set up using memory allocation. (N1)
	for (int c = 0; c < NUM_CELLS; c++) {
		atomsPerCell[c] = 0;
	}

	// Set atoms into layers (N2)
	for (int i = 0; i < NUM_ATOMS; i++) {
		// Randomly select a cell in atomsPerCell to put this atom in
		int cellIndex = rand() % NUM_CELLS;
		atomsPerCell[cellIndex] += 1;
		// Store the layer that this atom is at in the cell
		atomsPerLayer[i] = atomsPerCell[cellIndex];
	}

	// Find the maximum number of layers (i.e., max(k_c) for c in [0, NUM_CELLS))
	// (N3)
	int numLayers = maxVal(atomsPerCell, NUM_CELLS);


	// 3. Create an OpenCL context:

	Context = clCreateContext( NULL, 1, &Device, NULL, NULL, &status );
	if( status != CL_SUCCESS )
		fprintf( stderr, "clCreateContext failed\n" );


	// 4. Create an OpenCL command queue:

	CmdQueue = clCreateCommandQueue( Context, Device, 0, &status );
	if( status != CL_SUCCESS )
		fprintf( stderr, "clCreateCommandQueue failed\n" );


	// 5. Allocate the GPU device memory buffers for the A, B and C matrices:

	size_t laticeSize = NUM_ATOMS * NUM_ATOMS * sizeof(float);
	int mw = NUM_ATOMS;
	size_t mwSize = sizeof(mw);
	size_t cSize = NUM_ATOMS * NUM_ATOMS * sizeof(float);

	// Allocating device memory for the atomic latice matrix
	cl_mem deviceLatice = clCreateBuffer( Context, CL_MEM_READ_ONLY, laticeSize, NULL, &status );
	if( status != CL_SUCCESS )
		fprintf( stderr, "clCreateBuffer failed for dLatice (1)\n" );
	
	// Allocating device memory for the width of the matrices
	cl_mem dMW = clCreateBuffer( Context, CL_MEM_READ_ONLY, mwSize, NULL, &status );
	if( status != CL_SUCCESS )
		fprintf( stderr, "clCreateBuffer failed for dMW (1)\n" );
	
	// Allocating device memory for the C matrix
	cl_mem dC = clCreateBuffer( Context, CL_MEM_WRITE_ONLY, cSize, NULL, &status );
	if( status != CL_SUCCESS )
		fprintf( stderr, "clCreateBuffer failed for dC (1)\n" );


	// 6. Enqueue the 3 commands to write the data from the host buffers to the device buffers:

	// Enqueue the data from matrix A to the device.
	status = clEnqueueWriteBuffer( CmdQueue, deviceLatice, CL_FALSE, 0, laticeSize, hLatice, 0, NULL, NULL );
	if( status != CL_SUCCESS )
		fprintf( stderr, "clEnqueueWriteBuffer failed for deviceLatice martix (1)\n" );
	
	// Enqueue the matrix width data dMW the device.
	status = clEnqueueWriteBuffer( CmdQueue, dMW, CL_FALSE, 0, mwSize, &mw, 0, NULL, NULL );
	if( status != CL_SUCCESS )
		fprintf( stderr, "clEnqueueWriteBuffer failed for matrix width mw (1)\n" );	

	Wait( CmdQueue );


	// This code is for the Molecular Dynamics model
	// 7. Read the kernel code from the molecular_dynamics.cl file ...

	// Read molecular_dynamics.cl for 
	fseek( fp, 0, SEEK_END );
	size_t fileSize = ftell( fp );
	fseek( fp, 0, SEEK_SET );
	char *clProgramTextMD = new char[ fileSize+1 ];		// leave room for '\0'
	size_t n = fread( clProgramTextMD, 1, fileSize, fp );
	clProgramTextMD[fileSize] = '\0';
	fclose( fp );
	if( n != fileSize )
		fprintf( stderr, "Expected to read %d bytes from '%s' -- actually read %d.\n", fileSize, CL_FILE_NAME, n );

	// ... and create the kernel program:

	char *strings[NUM_MODELS];
	strings[0] = clProgramTextMD;	// Add the MD kernel to the list of string pointers for use with creating the program.
	Program = clCreateProgramWithSource( Context, NUM_MODELS, (const char **)strings, NULL, &status );	// IMPORTANT: Multiple kernel .cl files can be read in.
	if( status != CL_SUCCESS )
		fprintf( stderr, "clCreateProgramWithSource failed\n" );
	delete [ ] clProgramTextMD;


	// 8. Compile and link the kernel code:

	char *options = { (char *)"" };
	status = clBuildProgram( Program, 1, &Device, options, NULL, NULL );
	if( status != CL_SUCCESS )
	{
		size_t size;
		clGetProgramBuildInfo( Program, Device, CL_PROGRAM_BUILD_LOG, 0, NULL, &size );
		cl_char *log = new cl_char[ size ];
		clGetProgramBuildInfo( Program, Device, CL_PROGRAM_BUILD_LOG, size, log, NULL );
		fprintf( stderr, "clBuildProgram failed:\n%s\n", log );
		delete [ ] log;
	}


	// 9. Create the kernel object for MD:

	// Create a Kernel for the MD function.
	Kernel = clCreateKernel( Program, "MolecularDynamics", &status );
	if( status != CL_SUCCESS )
		fprintf( stderr, "clCreateKernel failed for MD\n" );


	// 10. setup the arguments to the kernel object:
	
	// Passing the dA matrix argument to the kernal object
	status = clSetKernelArg( Kernel, 0, sizeof(cl_mem), &deviceLatice);
	if ( status != CL_SUCCESS )
		fprintf( stderr, "clSetKernelArg failed for deviceLatice (%d)\n", status);

	// Passing the dMW matrix width argument to the kernal object
	status = clSetKernelArg( Kernel, 2, sizeof(cl_mem), &dMW);
	if ( status != CL_SUCCESS )
		fprintf( stderr, "clSetKernelArg failed for dMW (%d)\n", status);

	// Passing the dC matrix argument to the kernal object
	status = clSetKernelArg( Kernel, 3, sizeof(cl_mem), &dC);
	if ( status != CL_SUCCESS )
		fprintf( stderr, "clSetKernelArg failed for dC (%d)\n", status);
	

	// 11. enqueue the kernel object for execution:

	size_t globalWorkSize[3] = { NUM_ATOMS,      NUM_ATOMS,      1 };
	size_t localWorkSize[3]  = { LOCALSIZE, LOCALSIZE, 1 };

#ifndef CSV
	fprintf( stderr, "Molecular Dynamics Model\n");
	fprintf( stderr, "Number of Work Groups = %5d x %5d\n", NUM_ATOMS/LOCALSIZE, NUM_ATOMS/LOCALSIZE );
#endif

	Wait( CmdQueue );

	double time0 = omp_get_wtime( );

	status = clEnqueueNDRangeKernel( CmdQueue, Kernel, 1, NULL, globalWorkSize, localWorkSize, 0, NULL, NULL);
	if( status != CL_SUCCESS )
		fprintf( stderr, "clEnqueueNDRangeKernel failed: %d\n", status );

	Wait( CmdQueue );
	double time1 = omp_get_wtime( );


	// 12. read the results buffer back from the device to the host:

	status = clEnqueueReadBuffer( CmdQueue, dC, CL_FALSE, 0, cSize, hC, 0, NULL, NULL );
	if( status != CL_SUCCESS )
			fprintf( stderr, "clEnqueueReadBuffer failed\n" );

	Wait( CmdQueue );


#ifdef CSV
	fprintf( stderr, "%8d , %6d , %10.2lf, %12.2f\n",
		NUM_ATOMS*NUM_ATOMS, LOCALSIZE*LOCALSIZE, (double)NUM_ATOMS*(double)NUM_ATOMS*(double)NUM_ATOMS/(time1-time0)/1000000000., hC[NUM_ATOMS-1][NUM_ATOMS-1] );
#else
	fprintf( stderr, "Molecular Dynamics Results");		// For MatrixMult, dC[NUM_ATOMS-1][NUM_ATOMS-1] = 2.0
	fprintf( stderr, "Matrix Size: %6d x %6d , Work Elements: %4d x %4d , GigaMultsPerSecond: %10.2lf, dC[%6d][%6d] = %12.2f\n",
		NUM_ATOMS, NUM_ATOMS, LOCALSIZE, LOCALSIZE, (double)NUM_ATOMS*(double)NUM_ATOMS*(double)NUM_ATOMS/(time1-time0)/1000000000., NUM_ATOMS-1, NUM_ATOMS-1, hC[NUM_ATOMS-1][NUM_ATOMS-1] );
#endif
	fprintf(stderr, "\n");


	// 13. clean everything up:

	clReleaseKernel(        Kernel   );
	clReleaseProgram(       Program  );
	clReleaseCommandQueue(  CmdQueue );
	clReleaseMemObject(     deviceLatice  );
	clReleaseMemObject(     dMW  );
	clReleaseMemObject(     dC  );

	return 0;
}


// Return the maximul value of an array of integers using binary reduction
int maxVal(int* s, int n) {
	int i;

	// Initialize maximum element
	int maximum = s[0];

	// Iterate through all of the elements to locate any values that might be larger
	// than the current maximum.
	for (i = 1; i < n; i++)
		if (s[i] > maximum)
			maximum = s[i];

	return maximum;
}

// wait until all queued tasks have taken place:

void Wait( cl_command_queue queue )
{
      cl_event wait;
      cl_int      status;

      status = clEnqueueMarker( queue, &wait );
      if( status != CL_SUCCESS )
	      fprintf( stderr, "Wait: clEnqueueMarker failed\n" );

      status = clWaitForEvents( 1, &wait );
      if( status != CL_SUCCESS )
	      fprintf( stderr, "Wait: clWaitForEvents failed\n" );
}


// vendor ids:
#define ID_AMD		0x1002
#define ID_INTEL	0x8086
#define ID_NVIDIA	0x10de

void SelectOpenclDevice()
{
		// select which opencl device to use:
		// priority order:
		//	1. a gpu
		//	2. an nvidia or amd gpu
		//	3. an intel gpu
		//	4. an intel cpu

	int bestPlatform = -1;
	int bestDevice = -1;
	cl_device_type bestDeviceType;
	cl_uint bestDeviceVendor;
	cl_int status;		// returned status from opencl calls
				// test against CL_SUCCESS

	// find out how many platforms are attached here and get their ids:

	cl_uint numPlatforms;
	status = clGetPlatformIDs(0, NULL, &numPlatforms);
	if( status != CL_SUCCESS )
		fprintf(stderr, "clGetPlatformIDs failed (1)\n");

	cl_platform_id* platforms = new cl_platform_id[numPlatforms];
	status = clGetPlatformIDs(numPlatforms, platforms, NULL);
	if( status != CL_SUCCESS )
		fprintf(stderr, "clGetPlatformIDs failed (2)\n");

	for( int p = 0; p < (int)numPlatforms; p++ )
	{
		// find out how many devices are attached to each platform and get their ids:

		cl_uint numDevices;

		status = clGetDeviceIDs(platforms[p], CL_DEVICE_TYPE_ALL, 0, NULL, &numDevices);
		if( status != CL_SUCCESS )
			fprintf(stderr, "clGetDeviceIDs failed (2)\n");

		cl_device_id* devices = new cl_device_id[numDevices];
		status = clGetDeviceIDs(platforms[p], CL_DEVICE_TYPE_ALL, numDevices, devices, NULL);
		if( status != CL_SUCCESS )
			fprintf(stderr, "clGetDeviceIDs failed (2)\n");

		for( int d = 0; d < (int)numDevices; d++ )
		{
			cl_device_type type;
			cl_uint vendor;
			size_t sizes[3] = { 0, 0, 0 };

			clGetDeviceInfo(devices[d], CL_DEVICE_TYPE, sizeof(type), &type, NULL);

			clGetDeviceInfo(devices[d], CL_DEVICE_VENDOR_ID, sizeof(vendor), &vendor, NULL);

			// select:

			if( bestPlatform < 0 )		// not yet holding anything -- we'll accept anything
			{
				bestPlatform = p;
				bestDevice = d;
				Platform = platforms[bestPlatform];
				Device = devices[bestDevice];
				bestDeviceType = type;
				bestDeviceVendor = vendor;
			}
			else					// holding something already -- can we do better?
			{
				if( bestDeviceType == CL_DEVICE_TYPE_CPU )		// holding a cpu already -- switch to a gpu if possible
				{
					if( type == CL_DEVICE_TYPE_GPU )			// found a gpu
					{										// switch to the gpu we just found
						bestPlatform = p;
						bestDevice = d;
						Platform = platforms[bestPlatform];
						Device = devices[bestDevice];
						bestDeviceType = type;
						bestDeviceVendor = vendor;
					}
				}
				else										// holding a gpu -- is a better gpu available?
				{
					if( bestDeviceVendor == ID_INTEL )			// currently holding an intel gpu
					{										// we are assuming we just found a bigger, badder nvidia or amd gpu
						bestPlatform = p;
						bestDevice = d;
						Platform = platforms[bestPlatform];
						Device = devices[bestDevice];
						bestDeviceType = type;
						bestDeviceVendor = vendor;
					}
				}
			}
		}
		delete [ ] devices;
	}
	delete [ ] platforms;


	if( bestPlatform < 0 )
	{
		fprintf(stderr, "I found no OpenCL devices!\n");
		exit( 1 );
	}
	else
	{
#ifndef CSV
		fprintf(stderr, "I have selected Platform #%d, Device #%d: ", bestPlatform, bestDevice);
		fprintf(stderr, "Vendor = %s, Type = %s\n", Vendor(bestDeviceVendor), Type(bestDeviceType) );
#endif
	}
}

char * Vendor( cl_uint v )
{
	switch( v )
	{
		case ID_AMD:
			return (char *)"AMD";
		case ID_INTEL:
			return (char *)"Intel";
		case ID_NVIDIA:
			return (char *)"NVIDIA";
	}
	return (char *)"Unknown";
}

char * Type( cl_device_type t )
{
	switch( t )
	{
		case CL_DEVICE_TYPE_CPU:
			return (char *)"CL_DEVICE_TYPE_CPU";
		case CL_DEVICE_TYPE_GPU:
			return (char *)"CL_DEVICE_TYPE_GPU";
		case CL_DEVICE_TYPE_ACCELERATOR:
			return (char *)"CL_DEVICE_TYPE_ACCELERATOR";
	}
	return (char *)"Unknown";
}