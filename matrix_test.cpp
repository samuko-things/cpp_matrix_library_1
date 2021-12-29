#include <iostream>
#include <stdio.h>
#include "matrixlab/matrixlab.h"

using namespace std;


int main()
{
   // float B[3][3]{
   //    {3,2,5},
   //    {4,7,9},
   //    {1,8,6},
   // };

   // float C[3][3];

   float A[3][3]{
      {2,-1,3},
      {1,3,-1},
      {2,-2,5},
   };

   float vect_A[3] = {2,11,3};
   float result[3];

   

   mat.print(A);

   mat.inverse(A, A); // A = A_inverse
   mat.print(A);

   mat.round(A, A,2); // A = round(A) to 2dp
   mat.print(A);

   mat.dot(result, A, vect_A); // result = A * vect_A (matrix dot multiplication)
   vect.print(result);

   vect.round(result, result, 0); // result = round(result) to the nearest whole num
   vect.print(result);
   

   return 0;
}
















// #include <iostream>
// #include "matrix_library/matrixlab.h"

// #include <stdio.h>
// #include <errno.h>

// using namespace std;


// int main(){
//     // Matrix A;
//     // float a[2][4]{
//     //     {2,3,5,1},
//     //     {4,6,0,7},
//     // };
//     // A.copyArray(a);

//     // Matrix B;
//     // float b[4][1]{
//     //     {3},
//     //     {4},
//     //     {2},
//     //     {9},
//     // };
//     // B.copyArray(b);
//     Matrix A;
//     float a[3][3]{
//         {1,2,3},
//         {4,1,5},
//         {6,0,2},
//     };
//     A.copyArray(a);

//     Matrix B;
//     float b[3][3]{
//         {2,7,4},
//         {3,1,6},
//         {5,0,8},
//     };
//     B.copyArray(b);

//     Matrix C = mat.dotProd(A, A.inverse());
//     C = C.round(0);
//     C.print();


//     // float detA = A.det();
//     // cout << "det of A = " << detA << endl;
  
//     // A = A.inverse();
//     // A.print();

//     // float detB = B.det();
//     // cout << "det of B = " << detB << endl;

//     B = B.inverse();
//     B.print();
//     B = B.round(3);
//     B.print();

//     return 0;
// }






// int add (){
//    return 0;
// }

// template <typename T, typename... ArgTypes>
// int add (T t, ArgTypes... args)
//  { return t + add(args...); }



