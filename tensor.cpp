#include <iostream>
#include <string>
#include <random>
#include <math.h>
#include <fstream>

#include "dais_exc.h"
#include "tensor.h"

#define PI 3.141592654
#define FLT_MAX 3.402823466e+38F /* max value */
#define FLT_MIN 1.175494351e-38F /* min positive value */

using namespace std;


/**
 * Random Initialization
 * 
 * Perform a random initialization of the tensor
 * 
 * @param mean The mean
 * @param std  Standard deviation
 */
void Tensor::init_random(float mean, float std){
    if(data){

        std::default_random_engine generator;
        std::normal_distribution<float> distribution(mean,std);

        for(int i=0;i<r;i++){
            for(int j=0;j<c;j++){
                for(int k=0;k<d;k++){
                    this->operator()(i,j,k)= distribution(generator);
                }
            }
        }    

    }else{
        throw(tensor_not_initialized());
    }
}


ostream& operator<< (ostream& is, const Tensor & t)
{
    if(t.r==0 || t.c==0 || t.d==0) is<<"Il tensore non e' inizzializzato";

    for(int k = 0; k < t.d; k++)//scorro per la profondita'sore non e' inizzializzato";
    {
        is <<"LIVELLO " << k << endl;//aggingo un identificati che mi definische a che profontità sono
        for(int j = 0; j < t.c; j++)//scorro per la colonne
        {
            for(int i = 0; i < t.r; i++)//scorro per la righe
                is << t(i, j, k) << "\t";     //aggiungo il valore allo stream e lo tabulo per seguire le colonne

            is << endl;//vado a capo ad ogni colonna
        }
        is << endl;//vado a capo ad ogni livello di profondita e salto una riga
    }
}

/**
 * Operator overloading - 
 * 1
 * It performs the point-wise difference between a Tensor and a constant
 * 
 * result(i,j,k)=this(i,j,k)-rhs
 * 
 * @return returns a new Tensor containing the result of the operation
 */
Tensor Tensor::operator-(const float &rhs) const
{
    Tensor a{(*this)};

    if(r==0 || c==0 || d==0) throw tensor_not_initialized{};

    for(int i = 0; i<a.r; i++)                             //scorro per le righe
        for(int j = 0; j<a.c; j++)                        //scorro per le colonne
            for(int k = 0; k<a.d; k++)a(i, j, k) -= rhs; //scorro per la profontita'

    return a;
}

/**
 * Operator overloading +
 * 
 * It performs the point-wise sum between a Tensor and a constant
 * 
 * result(i,j,k)=this(i,j,k)+rhs
 * 
 * @return returns a new Tensor containing the result of the operation
 */
Tensor Tensor::operator+(const float &rhs) const
{
    Tensor a{(*this)};

    if(r==0 || c==0 || d==0) throw tensor_not_initialized{};//lancio errore di vettore non inizzializzato 

    for(int i = 0; i<a.r; i++)                             //scorro per le righe
        for(int j = 0; j<a.c; j++)                        //scorro per le colonne
            for(int k = 0; k<a.d; k++)a(i, j, k) += rhs; //scorro per la profontita'

    return a;
}

/**
 * Operator overloading *
 * 
 * It performs the point-wise product between a Tensor and a constant
 * 
 * result(i,j,k)=this(i,j,k)*rhs
 * 
 * @return returns a new Tensor containing the result of the operation
 */
Tensor Tensor::operator*(const float &rhs) const
{
    Tensor a{(*this)};
    if(r==0 || c==0 || d==0) throw tensor_not_initialized{};//lancio errore di vettore non inizzializzato 

    for(int i = 0; i<a.r; i++)                             //scorro per le righe
        for(int j = 0; j<a.c; j++)                        //scorro per le colonne
            for(int k = 0; k<a.d; k++)a(i, j, k) *= rhs; //scorro per la profontita

    return a;
}

/**
 * Operator overloading / between a Tensor and a constant
 * 
 * It performs the point-wise division between a Tensor and a constant
 * 
 * result(i,j,k)=this(i,j,k)/rhs
 * 
 * @return returns a new Tensor containing the result of the operation
 */
Tensor Tensor::operator/(const float &rhs) const
{
    Tensor a{(*this)};

    if(r==0 || c==0 || d==0) throw tensor_not_initialized{};//lancio errore di vettore non inizzializzato 
    else if(rhs == 0)        throw unknown_operation{};//lancio errore di operazione sconosciuta

    for(int i = 0; i<a.r; i++)                             //scorro per le righe
        for(int j = 0; j<a.c; j++)                        //scorro per le colonne
            for(int k = 0; k<a.d; k++)a(i, j, k) /= rhs; //scorro per la profontita'

    return a;
}

/**
 * Operator overloading = (assignment) 
 * 
 * Perform the assignment between this object and another
 * 
 * @return a reference to the receiver object
 */
Tensor & Tensor::operator=(const Tensor &other)
{
    int n = other.r * other.c * other.d;

    this->r = other.r;
    this->c = other.c;
    this->d = other.d;

    if(this->data != other.data){
        if(this->r < 0 && this->c < 0 && this->d < 0)
            throw unknown_operation();

        if(data)
            this->~Tensor();
        
        data = new float[n];

        for(int i = 0; i < n; i++){
            this->data[i] = other.data[i];
        }
    }
    return (*this);
}

/**
 * Tensor padding
 * 
 * Zero pad a tensor in height and width, the new tensor will have the following dimensions:
 * 
 * (rows+2*pad_h) x (cols+2*pad_w) x (depth) 
 * 
 * @param pad_h the height padding
 * @param pad_w the width padding
 * @return the padded tensor
 */
Tensor Tensor::padding(int pad_h, int pad_w) const
{
    if(r==0 || c==0 || d==0) throw tensor_not_initialized{};//lancio errore di vettore non inizzializzato 

    Tensor rit{r+pad_h*2, c+pad_w*2, d, 0}; //creo un nuovo oggetto di tipo tensore copiando quello passato

    for(int i = pad_h; i<rit.r-pad_h; i++)      //scorro per le righe
        for(int j = pad_w; j<rit.c-pad_w; j++) //scorro per le colonne
            for(int k = 0; k<rit.d; k++)      //scorro per la profondita' 
                rit(i, j, k) = (*this)(i-pad_h, j-pad_w, k); //sostituisco nel tensore di ritorno il valore del tensore su cui e' stata chiamata la funzione
    return rit;
}

/**
 * Subset a tensor
 * 
 * retuns a part of the tensor having the following indices:
 * row_start <= i < row_end  
 * col_start <= j < col_end 
 * depth_start <= k < depth_end
 * 
 * The right extrema is NOT included
 * 
 * @param row_start 
 * @param row_end 
 * @param col_start
 * @param col_end
 * @param depth_start
 * @param depth_end
 * @return the subset of the original tensor
 */
Tensor Tensor::subset(unsigned int row_start, unsigned int row_end, unsigned int col_start, unsigned int col_end, unsigned int depth_start, unsigned int depth_end) const
{
    Tensor rit(row_end-row_start, col_end-col_start, depth_end-depth_start, 0); //creo un nuovo oggetto di tipo sensore copiando quello passato

    if(r==0 || c==0 || d==0)               throw tensor_not_initialized{};//lancio errore di vettore non inizzializzato 
    else if(rit.c>c || rit.r>r || rit.d>d) throw dimension_mismatch{};   //lancio errore di dimensioni non corrispondenti 

    for(size_t i = row_start; i<row_end; i++)               //scorro per le righe
        for(size_t j = col_start; j<col_end; j++)          //scorro per le colonne
            for(size_t k = depth_start; k<depth_end; k++) //scorro per la profondita' 
                rit(i-row_start, j-col_start, k-depth_start) = (*this)(i, j, k); //sostituisco nel tensore di ritorno il valore del tensore su cui e' stata chiamata la funzione
    return rit;
}