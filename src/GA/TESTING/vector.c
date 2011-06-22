/*-----------------------------------------------------------------------
/
/   vector.c routine for handling vectors
/
/------------------------------------------------------------------------*/
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "vector.h"
#define PRINT_VEC 1
/*-----------------------------------------------------------
/
/   vec_diff - take the difference between two vectors
/
/------------------------------------------------------------*/
vector vec_diff(vector* a2, vector* a1)
{
  vector a3;

  a3.x = a2->x - a1->x;
  a3.y = a2->y - a1->y;
  a3.z = a2->z - a1->z;
  return a3;
}
/*-----------------------------------------------------------
/
/   vec_sum - take the sum of two vectors
/
/------------------------------------------------------------*/
vector vec_sum(vector* a2,vector* a1)
{
  vector a3;

  a3.x = a2->x + a1->x;
  a3.y = a2->y + a1->y;
  a3.z = a2->z + a1->z;
  return a3;
}
/*-----------------------------------------------------------
/
/   vec_mult - multiply a vector by a scalar
/
/------------------------------------------------------------*/
vector vec_mult(vector* a1, double scalar)
{
  vector a2;

  a2.x = a1->x * scalar;
  a2.y = a1->y * scalar;
  a2.z = a1->z * scalar;
  return a2;
}
/*-----------------------------------------------------------
/
/   vec_3cm - return a vector that is in the center of a triangle 
/
/------------------------------------------------------------*/
vector vec_3cm(vector* a1, vector* a2, vector* a3)
{
  vector acm;

  acm.x = (a1->x + a2->x + a3->x)/3.0;
  acm.y = (a1->y + a2->y + a3->y)/3.0;
  acm.z = (a1->z + a2->z + a3->z)/3.0;
  return acm;
}
/*-----------------------------------------------------------
/
/   dot_prod - calculate the dot product of two vectors
/
/------------------------------------------------------------*/
double dot_prod(vector* a1,vector* a2)
{
  double product;

  product = a2->x * a1->x + a2->y * a1->y +a2->z * a1->z;
  return product;
}
/*-----------------------------------------------------------
/
/   cross_prod - calculate the cross product of two vectors
/
/------------------------------------------------------------*/
vector cross_prod(vector* a1,vector* a2)
{
  vector a3;

  a3.x = a1->y * a2->z - a1->z * a2->y; 
  a3.y = a1->z * a2->x - a1->x * a2->z; 
  a3.z = a1->x * a2->y - a1->y * a2->x; 
  return a3;
}
/*-----------------------------------------------------------
/
/   vec_norm - normalize a vector
/
/------------------------------------------------------------*/
void vec_norm(vector* a)
{
  double length;
  length = sqrt(a->x*a->x + a->y*a->y + a->z*a->z); 
  a->x/=length;
  a->y/=length;
  a->z/=length;
}
/*-----------------------------------------------------------
/
/   vec_zero - zero a vector
/
/------------------------------------------------------------*/
void vec_zero(vector* a)
{
  a->x = 0.0;
  a->y = 0.0;
  a->z = 0.0;
}
/*-----------------------------------------------------------
/
/   vec_length - get the length of a vector
/
/------------------------------------------------------------*/
double vec_length(vector* a)
{
  double length;
  length = sqrt(a->x*a->x + a->y*a->y + a->z*a->z); 
  return length;
}
/*-----------------------------------------------------------
/
/   vec_dist - get the distance between two vectors
/
/------------------------------------------------------------*/
double vec_dist(vector* a1, vector* a2)
{
  vector a12;
  double length;

  a12 = vec_diff(a1,a2);
  length = vec_length(&a12);
  return length;
}
/*-----------------------------------------------------------
/
/   vec_rtp - get r,theta,phi of a vector
/
/------------------------------------------------------------*/
void vec_rtp(vector* a,double* rv,double* ctv,double* stv,double* cpv,double* spv)
{
  int sign;

  (*rv)  = vec_length(a);
  (*ctv) = a->z/ (*rv);
  (*stv) = sqrt(a->x*a->x + a->y*a->y) / (*rv);

  if(fabs(*stv) <= 1.0e-8) {
    (*cpv) = 1.0;
    (*spv) = 0.0;
  }
  else {
    (*cpv) = a->x / (*rv) / (*stv);
    (*spv) = a->y / (*rv) / (*stv);
  }
}
/*-----------------------------------------------------------
/
/   vec_print - print a vector
/
/------------------------------------------------------------*/
void vec_print(vector* a,char* string)
{
  if(PRINT_VEC) 
    fprintf(stderr,"%s\t[%12.6f] [%12.6f] [%12.6f] \n",string,a->x,a->y,a->z);
}

