#include <stdlib.h>
#include <math.h>
#include "vector.h"
#include "config.h"
#include <cuda_runtime.h>

//compute: Updates the positions and locations of the objects in the system based on gravity.
//Parameters: None
//Returns: None
//Side Effect: Modifies the hPos and hVel arrays with the new positions and accelerations after 1 INTERVAL

__global__ void compute_kernel(vector3 *values, double *hPos, double *hVel, double *mass){
	int i,j,k;

	//	vector3** accels=(vector3**)malloc(sizeof(vector3*)*NUMENTITIES);
	// vector3* values=(vector3*)malloc(sizeof(vector3)*NUMENTITIES*NUMENTITIES);
	// for (i=0;i<NUMENTITIES;i++) // what is this doing
	//	accels[i]=&values[i*NUMENTITIES];
	
	
	//first compute the pairwise accelerations.  Effect is on the first argument.

	i = blockIdx.x * blockDim.x + threadIdx.x;
	j = blockIdx.y * blockDim.y + threadIdx.y;

	// must stay in bounds
	if (i >= NUMENTITIES || j >= NUMENTITIES){
		return;
	}
	// planet goes not have affect on itself
	if (i==j) {
		FILL_VECTOR(values[i*NUMENTITIES + j],0,0,0);
	}
	else{
		vector3 distance;
		for (k=0;k<3;k++) distance[k]=hPos[i*3 + k]-hPos[j*3 + k];
		double magnitude_sq=distance[0]*distance[0]+distance[1]*distance[1]+distance[2]*distance[2];
		double magnitude=sqrt(magnitude_sq);
		double accelmag=-1*GRAV_CONSTANT*mass[j]/magnitude_sq;
		FILL_VECTOR(values[i*NUMENTITIES + j],accelmag*distance[0]/magnitude,accelmag*distance[1]/magnitude,accelmag*distance[2]/magnitude);
	}
}

__global__ void add_kernel(vector3 *values, double *hPos, double *hVel, double *mass){

	int i = blockIdx.x * blockDim.x + threadIdx.x;
	int j, k;

	//sum up the rows of our matrix to get effect on each entity, then update velocity and position.
	
	if (i > NUMENTITIES) return;

	vector3 accel_sum={0,0,0};
	for (j=0;j<NUMENTITIES;j++){
		for (k=0;k<3;k++)
			accel_sum[k]+=values[i*NUMENTITIES + j][k];
	}

	//compute the new velocity based on the acceleration and time interval
	//compute the new position based on the velocity and time interval
	for (k=0;k<3;k++){
		hVel[i*3 + k]+=accel_sum[k]*INTERVAL;
		hPos[i*3 + k]+=hVel[i*3 + k]*INTERVAL;
	}
		
}


void compute(){
	
	// creating variables (host)
	vector3 *all_values;
	double *d_hPos, *d_hVel, *d_mass;
	size_t all_values_size = sizeof(vector3) * NUMENTITIES * NUMENTITIES;
	size_t pos_vel_size = sizeof(double) * NUMENTITIES * 3;
	size_t mass_size = sizeof(double) * NUMENTITIES;
   
	// allocating space for variables to pass to function (device)
	cudaMalloc((void **)&all_values, all_values_size);
	cudaMalloc((void **)&d_hPos, pos_vel_size);
	cudaMalloc((void **)&d_hVel, pos_vel_size);
	cudaMalloc((void **)&d_mass, mass_size);

	// copy data from host to device
	cudaMemcpy(d_hPos, hPos, pos_vel_size, cudaMemcpyHostToDevice);
	cudaMemcpy(d_hVel, hVel, pos_vel_size, cudaMemcpyHostToDevice);
	cudaMemcpy(d_mass, mass, mass_size, cudaMemcpyHostToDevice);

	// kernel call for compute
	dim3 dimBlock(16, 16);
	dim3 dimGrid((NUMENTITIES + 15)/16, (NUMENTITIES + 15)/16);
	compute_kernel<<<dimGrid, dimBlock>>>(all_values, d_hPos, d_hVel, d_mass);
	cudaDeviceSynchronize();

	// kernel call for add
	add_kernel<<<NUMENTITIES, 1>>>(all_values, d_hPos, d_hVel, d_mass);
	cudaDeviceSynchronize();

	// what is the reduciton portion?? is it just the addition

	// copy updated data device back to host
	cudaMemcpy(hPos, d_hPos, pos_vel_size, cudaMemcpyDeviceToHost);
	cudaMemcpy(hVel, d_hVel, pos_vel_size, cudaMemcpyDeviceToHost);
	// don't need to copy mass back

	// free memory of variables (device)
	cudaFree(all_values);
	cudaFree(d_hPos);
	cudaFree(d_hVel);
	cudaFree(d_mass);
	
	//free(accels);

}
