/*----------------------------------------------------------------
/
/   vector.h
/
/----------------------------------------------------------------*/
typedef struct {
  double x;
  double y;
  double z;
} vector;

vector vec_diff(vector*,vector*);
vector vec_sum(vector*,vector*);
vector vec_mult(vector*,double);
double dot_prod(vector*,vector*);
vector cross_prod(vector*,vector*);
void vec_norm(vector*);
void vec_zero(vector*);
double vec_length(vector*);
void vec_rtp(vector*,double*,double*,double*,double*,double*);
void vec_print(vector*,char*);
vector vec_3cm(vector*, vector*, vector*);
double vec_dist(vector*,vector*);
