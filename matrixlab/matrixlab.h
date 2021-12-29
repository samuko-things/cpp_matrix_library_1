
#ifndef _MATRIXLAB_H_
#define _MATRIXLAB_H_
#endif

#include <iostream>
#include <math.h>
using namespace std;

typedef unsigned int SIZE;



/* the rounding to decimal place algorithm was gotten from geeks for geeks and edited by me*/
float round_dp(float val, int dp){
    if (val>=0){
        float rounded_val = (int)((val*pow(10,dp))+0.5);
        return (float)rounded_val/pow(10,dp);
    }
    else{
        float rounded_val = (int)((val*pow(10,dp))-0.5);
        return (float)rounded_val/pow(10,dp);
    }
}










class VectorOperations{
    public:
        // this function prints the vector
        template<SIZE N>
        void print(float (&vector1)[N]) {
            cout << N << "D Vector" << endl;
            for (SIZE c = 0; c < N; c += 1) {
                if (c == N-1) cout << vector1[c];
                else cout << vector1[c] << ", ";
            }
            cout << endl;
            cout << endl;
        }

        void printWithDirection(float (&vector1)[3]) {
            cout << 3 << "D Vector" << endl;
            cout << vector1[0] << "i +" << vector1[1] << "j +" << vector1[2] << "k" << endl;
            cout << endl;
        }

        void printWithDirection(float (&vector1)[2]) {
            cout << 2 << "D Vector" << endl;
            cout << vector1[0] << "i +" << vector1[1] << "j" << endl;
            cout << endl;
        }


        // this function copies a vector to another vector of the same size
        template<SIZE N>
        void copy(float (&destinationVector)[N], float (&sourceVector)[N]) {
            for (SIZE c = 0; c < N; c += 1) {
                destinationVector[c] = sourceVector[c];
            }
        }


        // this function clears a vector
        template<SIZE N>
        void clear(float (&vector1)[N]) {
            for (SIZE c = 0; c < N; c += 1) {
                vector1[c] = 0;
            }
        }


        // this function rounds a vector to the nearest decimal place
        template<SIZE N>
        void round(float (&rounded_vector)[N], float (&vector1)[N], int dp) {
            float buffer[N];
            for (int c = 0; c < N; c += 1) {
                buffer[c] = round_dp(vector1[c], dp);
            }
            clear(rounded_vector);
            copy(rounded_vector, buffer);
        }


        // this function adds two vectors together
        template<SIZE N>
        void add(float (&result)[N], float (&vector1)[N], float (&vector2)[N]) {
            float buffer[N];
            for (SIZE c = 0; c < N; c += 1) {
                buffer[c] = vector1[c] + vector2[c];
            }
            clear(result);
            copy(result, buffer);
        }


        // this function subtracts a vector from another
        template<SIZE N>
        void subtract(float (&result)[N], float (&vector1)[N], float (&vector2)[N]) {
            float buffer[N];
            for (SIZE c = 0; c < N; c += 1) {
                buffer[c] = vector1[c] - vector2[c];
            }
            clear(result);
            copy(result, buffer);
        }


        // this function multiply two vectors together element by element
        template<SIZE N>
        void multiply(float (&result)[N], float (&vector1)[N], float (&vector2)[N]) {
            float buffer[N];
            for (SIZE c = 0; c < N; c += 1) {
                buffer[c] = vector1[c] * vector2[c];
            }
            clear(result);
            copy(result, buffer);
        }


        // this function scale a vector
        template<SIZE N>
        void scale(float (&result)[N], float (&vector1)[N], float num) {
            float buffer[N];
            for (SIZE c = 0; c < N; c += 1) {
                buffer[c] = vector1[c] * num;
            }
            clear(result);
            copy(result, buffer);
        }

        template<SIZE N>
        void scaleDiv(float (&result)[N], float (&vector1)[N], float num) {
            float buffer[N];
            for (SIZE c = 0; c < N; c += 1) {
                buffer[c] = vector1[c] / num;
            }
            clear(result);
            copy(result, buffer);
        }

    
        // this function performs dot product on two vector
        template<SIZE N>
        float dot(float (&vector1)[N], float (&vector2)[N]) {
            float dot_prod;
            for (int c = 0; c < N; c += 1) {
                dot_prod += (vector1[c] * vector2[c]);
            }
            return dot_prod;
        }


        // this function performs cross product on two vector
        void cross(float (&result)[3], float (&vector1)[2], float (&vector2)[2]) {
            float buffer[3];

            buffer[0] = (vector1[1] * 0) - (0 * vector2[1]);
            buffer[1] = -1 * ( (vector1[0] * 0) - (0 * vector2[0]) );
            buffer[2] = (vector1[0] * vector2[1]) - (vector1[1] * vector2[0]);

            clear(result);
            copy(result, buffer);
        }

        void cross(float (&result)[3], float (&vector1)[3], float (&vector2)[3]) {
            float buffer[3];

            buffer[0] = (vector1[1] * vector2[2]) - (vector1[2] * vector2[1]);
            buffer[1] = -1 * ( (vector1[0] * vector2[2]) - (vector1[2] * vector2[0]) );
            buffer[2] = (vector1[0] * vector2[1]) - (vector1[1] * vector2[0]);

            clear(result);
            copy(result, buffer);
        }


	    // this function calculates the magnitude of a vector
        template<SIZE N>
        float magnitude(float (&vector1)[N]) {
            float mag;
            for (int c = 0; c < N; c += 1) {
                mag += pow(vector1[c],2);
            }
            mag = sqrt(mag);

            return mag;
        }
	

        // this function normalize a vector (it generates unit vector for a 3D or 2D vector)
        template<SIZE N>
        void normalize(float (&result)[N], float (&vector1)[N]) {
            float buffer[N];
            float mag = magnitude(vector1);

            for (int c = 0; c < N; c += 1) {
                buffer[N] = vector1[N] / mag;
            }

            clear(result);
            copy(result, buffer);
        }


        // this function calculates the Direction Cosine in radians
        void directionCosinesRad(float (&result)[3], float (&vector1)[3]) {
            float buffer[3];
            float unit_vector[3];
            normalize(unit_vector, vector1);

            buffer[0] = acos(unit_vector[0]);
            buffer[1] = acos(unit_vector[1]);
            buffer[2] = acos(unit_vector[2]);

            clear(result);
            copy(result, buffer);
        }


        // this function calculates the Direction Cosine in radians
        void directionCosinesDeg(float (&result)[3], float (&vector1)[3]) {
            float buffer[3];
            float unit_vector[3];
            normalize(unit_vector, vector1);

            buffer[0] = acos(unit_vector[0])*180/M_PI;
            buffer[1] = acos(unit_vector[1])*180/M_PI;
            buffer[2] = acos(unit_vector[2])*180/M_PI;

            clear(result);
            copy(result, buffer);
        }



        // this function calculates the angle between two vectors
        float angleBtwRad(float (&vector1)[3], float (&vector2)[3]) {
            float unit_vector1[3], unit_vector2[3];

            normalize(unit_vector1, vector1);
            normalize(unit_vector2, vector2);

            float angle = acos( (unit_vector1[0] * unit_vector2[0]) + (unit_vector1[1] * unit_vector2[1]) + (unit_vector1[2] * unit_vector2[2]) );
            return angle;
        }


        // this function calculates the angle between two vectors
        float angleBtwDeg(float (&vector1)[3], float (&vector2)[3]) {
            float unit_vector1[3], unit_vector2[3];

            normalize(unit_vector1, vector1);
            normalize(unit_vector2, vector2);

            float angle = acos( (unit_vector1[0] * unit_vector2[0]) + (unit_vector1[1] * unit_vector2[1]) + (unit_vector1[2] * unit_vector2[2]) ) * 180 / M_PI;
            return angle;
        }



    private:

        float round_dp(float val, int dp){
            if (val>=0){
                float rounded_val = (int)((val*pow(10,dp))+0.5);
                return (float)rounded_val/pow(10,dp);
            }
            else{
                float rounded_val = (int)((val*pow(10,dp))-0.5);
                return (float)rounded_val/pow(10,dp);
            }
        }

};

VectorOperations vect;











class MatrixOperations{
    public:

        // this function prints the matrix
        template<SIZE R, SIZE C>
        void print(float (&matrix)[R][C]) {
            cout << R << " by " << C << " Matrix" << endl;
            for (SIZE r = 0; r < R; r += 1) {
                for (SIZE c = 0; c < C; c += 1) {
                cout << matrix[r][c] << "\t";
                }
                cout << endl;
            }
            cout << endl;
        }

        // this function copies a matrix to another matrix of the same size
        template<SIZE R, SIZE C>
        void copy(float (&destinationMatrix)[R][C], float (&sourceMatrix)[R][C]) {
            for (SIZE r = 0; r < R; r += 1) {
                for (SIZE c = 0; c < C; c += 1) {
                    destinationMatrix[r][c] = sourceMatrix[r][c];
                }
            }
        }


        // this function clears a matrix
        template<SIZE R, SIZE C>
        void clear(float (&matrix)[R][C]) {
            for (SIZE r = 0; r < R; r += 1) {
                for (SIZE c = 0; c < C; c += 1) {
                    matrix[r][c] = 0;
                }
            }
        }

        // this function adds two matrices or vectors together
        template<SIZE R, SIZE C>
        void add(float (&result)[R][C], float (&matrix1)[R][C], float (&matrix2)[R][C]) {
            float buffer[R][C];
            for (SIZE r = 0; r < R; r += 1) {
                for (SIZE c = 0; c < C; c += 1) {
                    buffer[r][c] = matrix1[r][c] + matrix2[r][c];
                }
            }
            clear(result);
            copy(result, buffer);
        }


        // this function subtracts a matrix or vector from another
        template<SIZE R, SIZE C>
        void subtract(float (&result)[R][C], float (&matrix1)[R][C], float (&matrix2)[R][C]) {
            float buffer[R][C];
            for (SIZE r = 0; r < R; r += 1) {
                for (SIZE c = 0; c < C; c += 1) {
                    buffer[r][c] = matrix1[r][c] - matrix2[r][c];
                }
            }
            clear(result);
            copy(result, buffer);
        }


        // this function transpose a matrix or vector
        template<SIZE R, SIZE C>
        void transpose(float (&T_matrix)[C][R], float (&matrix)[R][C]) {
            float buffer[C][R];
            for (SIZE r = 0; r < C; r += 1) {
                for (SIZE c = 0; c < R; c += 1) {
                    buffer[r][c] = matrix[c][r];
                }
            }
            clear(T_matrix);
            copy(T_matrix, buffer);
        }


        // this function multiplies two matrices together
        template<SIZE R, SIZE C, SIZE N>
        void dot(float (&result)[R][C], float (&matrix1)[R][N], float (&matrix2)[N][C]) {
            float buffer[R][C];
            float val;
            float T_matrix2[C][N];
            transpose(T_matrix2, matrix2); // transpose matrix2 and store ans in T_matrix2

            for (SIZE row = 0; row < R; row += 1) {
                for (SIZE col = 0; col < C; col += 1) {
                    for (SIZE count = 0; count < N; count += 1) {
                        val += matrix1[row][count] * T_matrix2[col][count];
                    }
                    buffer[row][col] = val;
                    val = 0;
                }
            }

            clear(result);
            copy(result, buffer);
        }


        // this function multiplies a matrix and a vector, stores the ans in a vector
        template<SIZE N>
        void dot(float (&result)[N], float (&matrix)[N][N], float (&vector1)[N]) {
            float buffer[N];
            float val;

            for (SIZE row = 0; row < N; row += 1) {
                for (SIZE col = 0; col < 1; col += 1) {
                    for (SIZE count = 0; count < N; count += 1) {
                        val += matrix[row][count] * vector1[count];
                    }
                    buffer[row] = val;
                    val = 0;
                }
            }

            vect.clear(result);
            vect.copy(result, buffer);
        }



        // this function multiplies a matrix or vector by a scaler
        template<SIZE R, SIZE C>
        void scale(float (&result)[R][C], float (&matrix1)[R][C], float num) {
            float buffer[R][C];
            for (SIZE r = 0; r < R; r += 1) {
                for (SIZE c = 0; c < C; c += 1) {
                    buffer[r][c] = matrix1[r][c]*num;
                }
            }
            clear(result);
            copy(result, buffer);
        }


        // this function divides a matrix or vector by a scaler
        template<SIZE R, SIZE C>
        void scaleDiv(float (&result)[R][C], float (&matrix1)[R][C], float num) {
            float buffer[R][C];
            for (SIZE r = 0; r < R; r += 1) {
                for (SIZE c = 0; c < C; c += 1) {
                    buffer[r][c] = matrix1[r][c]/num;
                }
            }
            clear(result);
            copy(result, buffer);
        }



        // this function multiplies matrices together element by element
        template<SIZE R, SIZE C>
        void multiply(float (&result)[R][C], float (&matrix1)[R][C], float (&matrix2)[R][C]) {
            float buffer[R][C];
            for (SIZE r = 0; r < R; r += 1) {
                for (SIZE c = 0; c < C; c += 1) {
                    buffer[r][c] = matrix1[r][c] * matrix2[r][c];
                }
            }
            clear(result);
            copy(result, buffer);
        }


        // this function calculates the determinant of a given matrix
        float det(float (&matrix)[1][1]) {	  
            return matrix[0][0];
        }

        template<SIZE N>
        float det(float (&matrix)[N][N]) {
            float d=0;
            float minor_matrix[N-1][N-1];

            for (SIZE c=0; c<N; c+=1){
                minor(minor_matrix, matrix,0,c);
                if(c==0 || (c%2)==0) d+=(matrix[0][c]*det(minor_matrix));
                else d-=(matrix[0][c]*det(minor_matrix));
            }
            return d;
        }


        // this generates inverse of a given matrix
        template<SIZE N>
        void inverse(float (&inv_matrix)[N][N], float (&matrix)[N][N]) {
            float d = det(matrix);

            float buffer[N][N];
            minorMatrix(buffer, matrix);
            cofactorMatrix(buffer, buffer);
            adjointMatrix(buffer, buffer);

            scaleDiv(inv_matrix, buffer, d); 
        }


        // this rounds a given matrix to the nearest decimal place , dp=0 gives the nearest whole number
        template<SIZE N>
        void round(float (&rounded_matrix)[N][N], float (&matrix)[N][N], int dp) {
            float buffer[N][N];
            for (int r = 0; r < N; r += 1) {
                for (int c = 0; c < N; c += 1) {
                    buffer[r][c] = round_dp(matrix[r][c], dp);
                }
            }
            clear(rounded_matrix);
            copy(rounded_matrix, buffer);
        }


    private:

        float round_dp(float val, int dp){
            if (val>=0){
                float rounded_val = (int)((val*pow(10,dp))+0.5);
                return (float)rounded_val/pow(10,dp);
            }
            else{
                float rounded_val = (int)((val*pow(10,dp))-0.5);
                return (float)rounded_val/pow(10,dp);
            }
        }
        

        // this functon gets the minor matrix element for a giving matrix position
        template<SIZE N>
        void minor(float (&minor_matrix)[N-1][N-1], float (&major_matrix)[N][N], int pos_r, int pos_c) {
            SIZE minor_r = 0, minor_c = 0;

            float buffer[N-1][N-1];

            if(N>=2){
                for(SIZE r=0; r<N; r+=1){

                    for(SIZE c=0; c<N; c+=1){
                        if((r==pos_r) || (c==pos_c)) {
                            continue;
                        }
                        else {
                            buffer[minor_r][minor_c] = major_matrix[r][c];
                            minor_c+=1;
                            if(minor_c>=N-1) minor_c = 0;
                        }
                    }

                if(r==pos_r) {
                            continue;
                    }
                    else {
                        minor_r+=1;
                        if(minor_r>=N-1) minor_r = 0;
                    } 
                }
            }

            clear(minor_matrix);
            copy(minor_matrix, buffer);
            
        }


        // this generates a new matrices of the det of minors of a given matrix
        template<SIZE N>
        void minorMatrix(float (&minor_matrix)[N][N], float (&major_matrix)[N][N]) {
            float buffer[N][N];
            float Minor[N-1][N-1];

            for(SIZE r=0; r<N; r+=1){
                for(SIZE c=0; c<N; c+=1){
                    minor(Minor, major_matrix,r,c);	
                    buffer[r][c] = det(Minor);
                }
            }
    
            clear(minor_matrix);
            copy(minor_matrix, buffer);
        }


        // this generates cofactor of a given matrix
        template<SIZE N>
        void cofactorMatrix(float (&co_matrix)[N][N], float (&minor_matrix)[N][N]) {
            float buffer[N][N];
            
            for(SIZE r=0; r<N; r+=1){
                for(SIZE c=0; c<N; c+=1){
                    if((r+c)==0 || ((r+c)%2)==0) buffer[r][c] = minor_matrix[r][c];
                    else buffer[r][c] = -1*minor_matrix[r][c];
                }
            }
    
            clear(co_matrix);
            copy(co_matrix, buffer);
        }


        // this generates adjoint of a given cofactor matrix
        template<SIZE N>
        void adjointMatrix(float (&ad_matrix)[N][N], float (&co_matrix)[N][N]) {
            transpose(ad_matrix, co_matrix);
        }


};

MatrixOperations mat;




































































// namespace mat{

//     // this function prints the matrix
// 	template<SIZE R, SIZE C>
// 	void print(float (&matrix)[R][C]) {
//         cout << R << " by " << C << " Matrix" << endl;
//         for (SIZE r = 0; r < R; r += 1) {
//             for (SIZE c = 0; c < C; c += 1) {
//             cout << matrix[r][c] << "\t";
//             }
//             cout << endl;
//         }
//         cout << endl;
// 	}

// 	// this function copies a matrix to another matrix of the same size
// 	template<SIZE R, SIZE C>
// 	void copy(float (&destinationMatrix)[R][C], float (&sourceMatrix)[R][C]) {
// 	  for (SIZE r = 0; r < R; r += 1) {
// 	    for (SIZE c = 0; c < C; c += 1) {
// 	      destinationMatrix[r][c] = sourceMatrix[r][c];
// 	    }
// 	  }float round(float val, int dp){
//     if (val>=0){
//         float rounded_val = (int)((val*pow(10,dp))+0.5);
//         return (float)rounded_val/pow(10,dp);
//     }
//     else{
//         float rounded_val = (int)((val*pow(10,dp))-0.5);
//         return (float)rounded_val/pow(10,dp);
//     }
// }
// 	}


// 	// this function clears a matrix
// 	template<SIZE R, SIZE C>
// 	void clear(float (&matrix)[R][C]) {
// 	  for (SIZE r = 0; r < R; r += 1) {
// 	    for (SIZE c = 0; c < C; c += 1) {
// 	      matrix[r][c] = 0;
// 	    }
// 	  }
// 	}

// 	// this function adds two matrices or vectors together
// 	template<SIZE R, SIZE C>
// 	void add(float (&result)[R][C], float (&matrix1)[R][C], float (&matrix2)[R][C]) {
// 	  float buffer[R][C];
// 	  for (SIZE r = 0; r < R; r += 1) {
// 	    for (SIZE c = 0; c < C; c += 1) {
// 	      buffer[r][c] = matrix1[r][c] + matrix2[r][c];
// 	    }
// 	  }

// 	  clear(result);
// 	  copy(result, buffer);
// 	}


// 	// this function subtracts a matrix or vector from another
// 	template<SIZE R, SIZE C>
// 	void subtract(float (&result)[R][C], float (&matrix1)[R][C], float (&matrix2)[R][C]) {
// 	  float buffer[R][C];
// 	  for (SIZE r = 0; r < R; r += 1) {
// 	    for (SIZE c = 0; c < C; c += 1) {
// 	      buffer[r][c] = matrix1[r][c] - matrix2[r][c];
// 	    }
// 	  }

// 	  clear(result);
// 	  copy(result, buffer);
// 	}


// 	// this function transpose a matrix or vector
// 	template<SIZE R, SIZE C>
// 	void transpose(float (&T_matrix)[C][R], float (&matrix)[R][C]) {
// 	  float buffer[C][R];
// 	  for (SIZE r = 0; r < C; r += 1) {
// 	    for (SIZE c = 0; c < R; c += 1) {
// 	      buffer[r][c] = matrix[c][r];
// 	    }
// 	  }

// 	  clear(T_matrix);
// 	  copy(T_matrix, buffer);
// 	}


// 	// this function multiplies two matrices together
// 	template<SIZE R, SIZE C, SIZE N>
// 	void dot(float (&result)[R][C], float (&matrix1)[R][N], float (&matrix2)[N][C]) {
// 	  float buffer[R][C];
// 	  float val;

// 	  float T_matrix2[C][N];
// 	  transpose(T_matrix2, matrix2); // transpose matrix2 and store ans in T_matrix2

// 	  for (SIZE row = 0; row < R; row += 1) {
// 	    for (SIZE col = 0; col < C; col += 1) {
// 	      for (SIZE count = 0; count < N; count += 1) {
// 	        val += matrix1[row][count] * T_matrix2[col][count];
// 	      }
// 	      buffer[row][col] = val;
// 	      val = 0;
// 	    }
// 	  }

// 	  clear(result);
// 	  copy(result, buffer);
// 	}



// 	// this function multiplies a matrix or vector by a scaler
// 	template<SIZE R, SIZE C>
// 	void scale(float (&result)[R][C], float (&matrix1)[R][C], float num) {
// 	  float buffer[R][C];
// 	  for (SIZE r = 0; r < R; r += 1) {
// 	    for (SIZE c = 0; c < C; c += 1) {
// 	      buffer[r][c] = matrix1[r][c]*num;
// 	    }
// 	  }

// 	  clear(result);
// 	  copy(result, buffer);
// 	}

//     // this function divides a matrix or vector by a scaler
//     template<SIZE R, SIZE C>
// 	void scaleDiv(float (&result)[R][C], float (&matrix1)[R][C], float num) {
// 	  float buffer[R][C];
// 	  for (SIZE r = 0; r < R; r += 1) {
// 	    for (SIZE c = 0; c < C; c += 1) {
// 	      buffer[r][c] = matrix1[r][c]/num;
// 	    }
// 	  }

// 	  clear(result);
// 	  copy(result, buffer);
// 	}



// 	// this function multiplies matrices together element by element
// 	template<SIZE R, SIZE C>
// 	void multiply(float (&result)[R][C], float (&matrix1)[R][C], float (&matrix2)[R][C]) {
// 	  float buffer[R][C];
// 	  for (SIZE r = 0; r < R; r += 1) {
// 	    for (SIZE c = 0; c < C; c += 1) {
// 	      buffer[r][c] = matrix1[r][c] * matrix2[r][c];
// 	    }
// 	  }
// 	  clear(result);
// 	  copy(result, buffer);
// 	}


//     // this functon gets the minor matrix element for a giving matrix position
//     template<SIZE N>
//     void minor(float (&minor_matrix)[N-1][N-1], float (&major_matrix)[N][N], int pos_r, int pos_c) {
//         SIZE minor_r = 0, minor_c = 0;

//         float buffer[N-1][N-1];

//         if(N>=2){
//             for(SIZE r=0; r<N; r+=1){

//                 for(SIZE c=0; c<N; c+=1){
//                     if((r==pos_r) || (c==pos_c)) {
//                         continue;
//                     }
//                     else {
//                         buffer[minor_r][minor_c] = major_matrix[r][c];
//                         minor_c+=1;
//                         if(minor_c>=N-1) minor_c = 0;
//                     }
//                 }

//             if(r==pos_r) {
//                         continue;
//                 }
//                 else {
//                     minor_r+=1;
//                     if(minor_r>=N-1) minor_r = 0;
//                 } 
//             }
//         }

//         clear(minor_matrix);
//         copy(minor_matrix, buffer);
	    
//     }



// 	// this function calculates the determinant of a given matrix
// 	float det(float (&matrix)[1][1]) {	  
// 	   return matrix[0][0];
// 	}

//     template<SIZE N>
//     float det(float (&matrix)[N][N]) {
//         float d=0;
//         float minor_matrix[N-1][N-1];

//         for (SIZE c=0; c<N; c+=1){
//             minor(minor_matrix, matrix,0,c);
//             if(c==0 || (c%2)==0) d+=(matrix[0][c]*det(minor_matrix));
//             else d-=(matrix[0][c]*det(minor_matrix));
//         }
//         return d;
//     }


//     // this generates a new matrices of the det of minors of a given matrix
//     template<SIZE N>
//     void minorMatrix(float (&minor_matrix)[N][N], float (&major_matrix)[N][N]) {
//         float buffer[N][N];
//         float Minor[N-1][N-1];

//         for(SIZE r=0; r<N; r+=1){
//             for(SIZE c=0; c<N; c+=1){
//                 minor(Minor, major_matrix,r,c);	
//                 buffer[r][c] = det(Minor);
//             }
//         }
   
//         clear(minor_matrix);
//         copy(minor_matrix, buffer);
//     }


//     // this generates cofactor of a given matrix
//     template<SIZE N>
//     void cofactorMatrix(float (&co_matrix)[N][N], float (&minor_matrix)[N][N]) {
//         float buffer[N][N];
        
//         for(SIZE r=0; r<N; r+=1){
//             for(SIZE c=0; c<N; c+=1){
//                 if((r+c)==0 || ((r+c)%2)==0) buffer[r][c] = minor_matrix[r][c];
//                 else buffer[r][c] = -1*minor_matrix[r][c];
//             }
//         }
   
//         clear(co_matrix);
//         copy(co_matrix, buffer);
//     }


//     // this generates adjoint of a given cofactor matrix
//     template<SIZE N>
//     void adjointMatrix(float (&ad_matrix)[N][N], float (&co_matrix)[N][N]) {
//         transpose(ad_matrix, co_matrix);
//     }

//     // this generates inverse of a given matrix
//     template<SIZE N>
//     void inverse(float (&inv_matrix)[N][N], float (&matrix)[N][N]) {
//         float d = det(matrix);

//         float buffer[N][N];
//         minorMatrix(buffer, matrix);
//         cofactorMatrix(buffer, buffer);
//         adjointMatrix(buffer, buffer);

//         scaleDiv(inv_matrix, buffer, d); 
//     }



//     // this rounds a given matrix to the nearest decimal place , dp=0 gives the nearest whole number
//     float round_dp(float val, int dp){
//         if (val>=0){
//             float rounded_val = (int)((val*pow(10,dp))+0.5);
//             return (float)rounded_val/pow(10,dp);
//         }
//         else{
//             float rounded_val = (int)((val*pow(10,dp))-0.5);
//             return (float)rounded_val/pow(10,dp);
//         }
//     }

//     template<SIZE N>
//     void round(float (&rounded_matrix)[N][N], float (&matrix)[N][N], int dp) {
//         float buffer[N][N];
//         for (int r = 0; r < N; r += 1) {
//             for (int c = 0; c < N; c += 1) {
//                 buffer[r][c] = round_dp(matrix[r][c], dp);
//             }
//         }
//         clear(rounded_matrix);
//         copy(rounded_matrix, buffer);
//     }

// }

























































































// #ifndef MATRIXLAB_h
// #define MATRIXLAB_h
// #endif

// #include <iostream>
// #include <math.h>
// using namespace std;


// const int R_SIZE = 5;
// const int C_SIZE = 5;


// class Matrix{
//     public:
//         int rowSize(){
//             return row;
//         }

//         int colSize(){
//             return col;
//         }

//         void writeSize(int r, int c){
//             __CLEAR();
//             row = r;
//             col = c;
//         }

//         void writeElement(int r_pos, int c_pos, float val){
//             if(r_pos<row && c_pos<col){
//                 data[r_pos][c_pos] = val;
//             }
//             else{
//                 cout << "no element exixst at this location" << endl;
//             }
//         }

//         float readElement(int r_pos, int c_pos){
//             if(r_pos<row && c_pos<col){
//                 return data[r_pos][c_pos];
//             }
//             else{
//                 cout << "no element exixst at this location" << endl;
//                 return 0;
//             }
//         }

//         Matrix round(int dp){
//             Matrix rounded_matrix;
//             rounded_matrix.row = row;
//             rounded_matrix.col = col;
//             for (int r = 0; r < row; r += 1) {
//                     for (int c = 0; c < col; c += 1) {
//                         rounded_matrix.data[r][c] = round_dp(data[r][c], dp);
//                     }
//             }
//             return rounded_matrix;
//         }

        
//         // copy c/c++ array into the matrix object
//         template<size_t R, size_t C>
//         void copyArray(float (&sourceArray)[R][C]) {
//             __CLEAR();
//             row = R; col = C;
//             for (int r = 0; r < R; r += 1) {
//                     for (int c = 0; c < C; c += 1) {
//                         data[r][c] = sourceArray[r][c];
//                     }
//             } 
//         }

//         // print the matrix
//         void print() {
//             cout << row << " by " << col << " Matrix" << endl;
//             for (int r = 0; r < row; r += 1) {
//                 for (int c = 0; c < col; c += 1) {
//                     cout << data[r][c] << '\t';
//                 }
//                 cout << endl;
//             }
//             cout << endl;
//         }

//         // transpose the matrix
//         Matrix Transpose() {     
//             Matrix transposed_matrix;         
//             transposed_matrix.row = col; transposed_matrix.col = row;
//             for (int r = 0; r < col; r += 1) {
//                 for (int c = 0; c < row; c += 1) {
//                     transposed_matrix.data[r][c] = data[c][r];
//                 }
//             }
//             return transposed_matrix;
//         }

//         Matrix scale(float scalar) {
//             Matrix scaled_matrix;
//             scaled_matrix.row = row; scaled_matrix.col = col;
//             for (int r = 0; r < row; r += 1) {
//                 for (int c = 0; c < col; c += 1) {
//                     scaled_matrix.data[r][c] = data[r][c] * scalar;
//                 }
//             }

//             return scaled_matrix;
//         }

//         Matrix scaleDiv(float scalar) {
//             Matrix scaled_matrix;
//             scaled_matrix.row = row; scaled_matrix.col = col;
//             for (int r = 0; r < row; r += 1) {
//                 for (int c = 0; c < col; c += 1) {
//                     scaled_matrix.data[r][c] = data[r][c] / scalar;
//                 }
//             }

//             return scaled_matrix;
//         }

        

//         float det() {
//             if(row == col){
//                 if (row ==1){
//                     det1();
//                 }
            
//                 else {
//                     float d=0;
//                     Matrix minor;
//                     for (int c=0; c<col; c+=1){
//                         minor = minorElement(0,c);
//                         if(c==0 || (c%2)==0) d+=(data[0][c]*minor.det());
//                         else d-=(data[0][c]*minor.det());
//                     }
//                     return d;
//                 }

            
//             }

//             else {
//                 cout << "invalid matrix determinant operation" << endl;
//                 return 0;
//             }   
//         }


//         Matrix inverse(){
//             Matrix inverse_matrix;
//             if(row == col){
//                 float d = det();
//                 inverse_matrix = minorMatrix();
//                 inverse_matrix = inverse_matrix.cofactorMatrix();
//                 inverse_matrix = inverse_matrix.adjointMatrix();
//                 inverse_matrix = inverse_matrix.scaleDiv(d);
//             }
//             else{
//                 cout << "wrong matrix size or dimension" << endl;
//             }

//            return inverse_matrix; 
//         }

// private:
//         int row, col;
//         float data[R_SIZE][C_SIZE];

//         float round_dp(float val, int dp){
//             float rounded_val = (int)((val*pow(10,dp))+0.5);
//             return (float)rounded_val/pow(10,dp);
//         }

//         // clear the whole matrix object
//         void __CLEAR() {
//             row = 0;col = 0;
//             for (int r = 0; r < R_SIZE; r += 1) {
//                 for (int c = 0; c < C_SIZE; c += 1) {
//                     data[r][c] = 0;
//                 }
//             }
//         }
        
//         float det1() {
//             if(row == col){
//                 float d=0;
//                 d = data[0][0];
//                 return d;
//             }  
//         }


//         Matrix minorElement(int pos_r, int pos_c) {
//             Matrix minor;
//             minor.row = row-1;
//             minor.col = col-1;

//             int minor_r = 0, minor_c = 0;

//             for(int r=0; r<row; r+=1){

//                 for(int c=0; c<col; c+=1){
//                     if((r==pos_r) || (c==pos_c)) {
//                         continue;
//                     }
//                     else {
//                         minor.data[minor_r][minor_c] = data[r][c];
//                         minor_c+=1;
//                         if(minor_c>=minor.col) minor_c = 0;
//                     }
//                 }

//                 if(r==pos_r) {
//                         continue;
//                 }
//                 else {
//                     minor_r+=1;
//                     if(minor_r>=minor.row) minor_r = 0;
//                 } 
//             }

//            return minor;  
//         }



//         Matrix minorMatrix(){
//             Matrix minor_element;
//             Matrix minor_matrix;
//             minor_matrix.row = row;
//             minor_matrix.col = col;

//             int minor_r = 0, minor_c = 0;

//             for(int r=0; r<row; r+=1){
//                 for(int c=0; c<col; c+=1){
//                     minor_element = minorElement(r,c);	
//                     minor_matrix.data[r][c] = minor_element.det();
//                 }
// 		    }

//            return minor_matrix; 
//         }


//         Matrix cofactorMatrix(){
//             Matrix cofactor_matrix;
//             cofactor_matrix.row = row;
//             cofactor_matrix.col = col;

//             for(int r=0; r<row; r+=1){
//                 for(int c=0; c<col; c+=1){
//                     if((r+c)==0 || ((r+c)%2)==0) cofactor_matrix.data[r][c] = data[r][c];
//                     else cofactor_matrix.data[r][c] = -1*data[r][c];
//                 }
//             }
//            return cofactor_matrix; 
//         }


//         Matrix adjointMatrix(){
//             Matrix adjoint_matrix;
//             adjoint_matrix = Transpose();
//            return adjoint_matrix; 
//         }
    
// };




// //////////////////////////// vector class ////////////////////////////////////
// class Vector{
//     public:
// };
// ///////////////////////////////////////////////////////////////////////////////





// ////////////////////////  matrices and vector operations ///////////////////////////////////////
// class MatrixOperations{
//     public:
//         Matrix add(Matrix matrix1, Matrix matrix2) {
//             Matrix add_matrix;
//             if((matrix1.rowSize() == matrix2.rowSize()) && (matrix1.colSize() == matrix2.colSize())){
//                 int new_row = matrix1.rowSize();
// 		        int new_col = matrix1.colSize();
//                 add_matrix.writeSize(new_row, new_col);

//                 float val;

//                 for (int r = 0; r < new_row; r += 1) {
//                     for (int c = 0; c < new_col; c += 1) {
//                         val = matrix1.readElement(r,c) + matrix2.readElement(r,c);
//                         add_matrix.writeElement(r,c,val);
//                         val = 0;
//                     }
//                 }
//                 return add_matrix;
//             }
//             else{
//                 cout << "ERROR: cannot add matrices, cross check their sizes" << endl;
//                 return add_matrix;
//             }
//         }



//         Matrix subtract(Matrix matrix1, Matrix matrix2) {
//             Matrix sub_matrix;
//             if((matrix1.rowSize() == matrix2.rowSize()) && (matrix1.colSize() == matrix2.colSize())){
//                 int new_row = matrix1.rowSize();
// 		        int new_col = matrix1.colSize();
//                 sub_matrix.writeSize(new_row, new_col);

//                 float val;

//                 for (int r = 0; r < new_row; r += 1) {
//                     for (int c = 0; c < new_col; c += 1) {
//                         val = matrix1.readElement(r,c) - matrix2.readElement(r,c);
//                         sub_matrix.writeElement(r,c,val);
//                         val = 0;
//                     }
//                 }
//                 return sub_matrix;
//             }
//             else{
//                 cout << "ERROR: cannot substract matrices, cross check their sizes" << endl;
//                 return sub_matrix;
//             }
//         }

        
//         Matrix dotProd(Matrix matrix1, Matrix matrix2) {
//             Matrix prod_matrix;
//             if(matrix1.colSize() == matrix2.rowSize()){
//                 int new_row = matrix1.rowSize();
// 		        int new_col = matrix2.colSize();
//                 prod_matrix.writeSize(new_row, new_col);

//                 Matrix T_matrix2;
//                 T_matrix2 = matrix2.Transpose();
//                 float val = 0;
//                 int loop = T_matrix2.colSize();

//                 for (int r=0; r<new_row; r+=1) {
//                     for (int c=0; c<new_col; c+=1) {
//                         for (int count = 0; count < loop; count += 1) {
//                             val += matrix1.readElement(r,count) * T_matrix2.readElement(c,count);
//                         }
//                         prod_matrix.writeElement(r,c,val);
//                         val = 0;
//                     }
//                 }
//                 return prod_matrix;
//             }
//             else{
//                 cout << "cannot multiply matrices, cross check their sizes" << endl;
//                 return prod_matrix;
//             }
//         }	
        
// };

// MatrixOperations mat;
// //////////////////////////////////////////////////////////////////////////////////////////////////////