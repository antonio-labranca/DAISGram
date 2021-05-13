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
 * Class constructor
 *
 * Parameter-less class constructor
 */
Tensor::Tensor() {}

/**
 * Class constructor
 *
 * Creates a new tensor of size r*c*d initialized at value v
 *
 * @param r
 * @param c
 * @param d
 * @param v
 * @return new Tensor
 */
Tensor::Tensor(int r, int c, int d, float v) {
	if (r < 0 || c < 0 || d < 0)
		throw unknown_operation();

	//Inizializzo righe, colonne e profondita'
	this->r = r;
	this->c = c;
	this->d = d;

	if (this->r != 0 && this->c != 0 && this->d != 0) {
		data = new float[r * c * d]; //Creo un vettore in memoria dinamica di dimensione righe x colonne x profondita'

		for (int i = 0; i < r; i++)
			for (int j = 0; j < c; j++)
				for (int k = 0; k < d; k++)
					(*this)(i, j, k) = v; //Inizializzo il tensore al valore v
	}
}

/**
 * Class distructor
 *
 * Cleanup the data when deallocated
 */
Tensor::~Tensor() {
	//Controllo se esiste il vettore data
	if (data) {
		delete[] data;
		data = nullptr;
	}
}

/**
 * Operator overloading ()
 *
 * if indexes are out of bound throw index_out_of_bound() exception
 *
 * @return the value at location [i][j][k]
 */
float Tensor::operator()(int i, int j, int k) const {
	//Controllo se le dimensioni del tensore sono in range [0, r] x [0, c] x [0, d]
	if (i >= 0 && i <= this->r && j >= 0 && j <= this->c && k >= 0 && k <= this->d)
		return data[(i * this->c + j) * this->d + k];
	else
		throw index_out_of_bound();
}

/**
 * Operator overloading ()
 *
 * Return the pointer to the location [i][j][k] such that the operator (i,j,k) can be used to
 * modify tensor data.
 *
 * If indexes are out of bound throw index_out_of_bound() exception
 *
 * @return the pointer to the location [i][j][k]
 */
float& Tensor::operator()(int i, int j, int k) {
	//Controllo se le dimensioni del tensore sono in range [0, r] x [0, c] x [0, d]
	if (i >= 0 && i <= this->r && j >= 0 && j <= this->c && k >= 0 && k <= this->d)
		return data[(i * this->c + j) * this->d + k];
	else
		throw index_out_of_bound();
}

/**
 * Copy constructor
 *
 * This constructor copies the data from another Tensor
 *
 * @return the new Tensor
 */
Tensor::Tensor(const Tensor& that) {
	//Inizializzo righe, colonne e profondita'
	this->r = that.r;
	this->c = that.c;
	this->d = that.d;

	//Controllo che il tensore that sia inizializzato
	if (that.r != 0 && that.c != 0 && that.d != 0) {
		//Creo un vettore di dimensioni righe x colonne x profondita'
		data = new float[r * c * d];

		for (int i = 0; i < r; i++)
			for (int j = 0; j < c; j++)
				for (int k = 0; k < d; k++)
					(*this)(i, j, k) = that(i, j, k);
	}
}

/**
 * Operator overloading -
 *
 * It performs the point-wise difference between two Tensors.
 *
 * result(i,j,k)=this(i,j,k)-rhs(i,j,k)
 *
 * The two tensors must have the same size otherwise throw a dimension_mismatch()
 *
 * @return returns a new Tensor containing the result of the operation
 */
Tensor Tensor::operator-(const Tensor& rhs)const {
	//controllo dimensioni
	if (this->r == rhs.r && this->c == rhs.c && this->d == rhs.d) {
		//creo nuovo tensore
		Tensor t(r, c, d);
		for (int i = 0; i < r; i++) {
			for (int j = 0; j < c; j++) {
				for (int k = 0; k < d; k++) {
					//svolgo operazione sottrazione
					t(i, j, k) = (*this)(i, j, k) - rhs(i, j, k);
				}
			}
		}
		return t;
	}
	else
		throw dimension_mismatch();
}

/**
 * Operator overloading +
 *
 * It performs the point-wise sum between two Tensors.
 *
 * result(i,j,k)=this(i,j,k)+rhs(i,j,k)
 *
 * The two tensors must have the same size otherwise throw a dimension_mismatch()
 *
 * @return returns a new Tensor containing the result of the operation
 */
Tensor Tensor::operator+(const Tensor& rhs)const {
	//controllo dimensioni
	if (this->r == rhs.r && this->c == rhs.c && this->d == rhs.d) {
		//creo nuovo tensore
		Tensor t(r, c, d);
		for (int i = 0; i < r; i++) {
			for (int j = 0; j < c; j++) {
				for (int k = 0; k < d; k++) {
					//svolgo operazione addizione
					t(i, j, k) = (*this)(i, j, k) + rhs(i, j, k);
				}
			}
		}
		return t;
	}
	else
		throw dimension_mismatch();
}

/**
 * Operator overloading *
 *
 * It performs the point-wise product between two Tensors.
 *
 * result(i,j,k)=this(i,j,k)*rhs(i,j,k)
 *
 * The two tensors must have the same size otherwise throw a dimension_mismatch()
 *
 * @return returns a new Tensor containing the result of the operation
 */
Tensor Tensor::operator*(const Tensor& rhs)const {
	//controllo dimensioni
	if (this->r == rhs.r && this->c == rhs.c && this->d == rhs.d) {
		//creo nuovo tensore
		Tensor t(r, c, d);
		for (int i = 0; i < r; i++) {
			for (int j = 0; j < c; j++) {
				for (int k = 0; k < d; k++) {
					//svolgo operazione moltiplicazione
					t(i, j, k) = (*this)(i, j, k) * rhs(i, j, k);
				}
			}
		}
		return t;
	}
	else
		throw dimension_mismatch();
}

/**
 * Operator overloading /
 *
 * It performs the point-wise division between two Tensors.
 *
 * result(i,j,k)=this(i,j,k)/rhs(i,j,k)
 *
 * The two tensors must have the same size otherwise throw a dimension_mismatch()
 *
 * @return returns a new Tensor containing the result of the operation
 */
Tensor Tensor::operator/(const Tensor& rhs) const {
	//controllo dimensioni
	if (this->r == rhs.r && this->c == rhs.c && this->d == rhs.d) {
		//creo nuovo tensore
		Tensor t(r, c, d);
		for (int i = 0; i < r; i++) {
			for (int j = 0; j < c; j++) {
				for (int k = 0; k < d; k++) {
					if (rhs(i, j, k) == 0)
						throw unknown_operation();
					// se non nullo svolgo  divisione
					t(i, j, k) = (*this)(i, j, k) / rhs(i, j, k);
				}
			}
		}
		return t;
	}
	else
		throw dimension_mismatch();
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
Tensor Tensor::operator-(const float& rhs) const
{
	Tensor a{ (*this) };

	if (r == 0 || c == 0 || d == 0) throw tensor_not_initialized{};

	for (int i = 0; i < a.r; i++)                             //scorro per le righe
		for (int j = 0; j < a.c; j++)                        //scorro per le colonne
			for (int k = 0; k < a.d; k++)a(i, j, k) -= rhs; //scorro per la profontita'

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
Tensor Tensor::operator+(const float& rhs) const
{
	Tensor a{ (*this) };

	if (r == 0 || c == 0 || d == 0) throw tensor_not_initialized{};//lancio errore di vettore non inizializzato 

	for (int i = 0; i < a.r; i++)                             //scorro per le righe
		for (int j = 0; j < a.c; j++)                        //scorro per le colonne
			for (int k = 0; k < a.d; k++)a(i, j, k) += rhs; //scorro per la profontita'

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
Tensor Tensor::operator*(const float& rhs) const
{
	Tensor a{ (*this) };
	if (r == 0 || c == 0 || d == 0) throw tensor_not_initialized{};//lancio errore di vettore non inizializzato 

	for (int i = 0; i < a.r; i++)                             //scorro per le righe
		for (int j = 0; j < a.c; j++)                        //scorro per le colonne
			for (int k = 0; k < a.d; k++)a(i, j, k) *= rhs; //scorro per la profontita'

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
Tensor Tensor::operator/(const float& rhs) const
{
	Tensor a{ (*this) };

	if (r == 0 || c == 0 || d == 0) throw tensor_not_initialized{};//lancio errore di vettore non inizializzato 
	else if (rhs == 0)        throw unknown_operation{};//lancio errore di operazione sconosciuta

	for (int i = 0; i < a.r; i++)                             //scorro per le righe
		for (int j = 0; j < a.c; j++)                        //scorro per le colonne
			for (int k = 0; k < a.d; k++)a(i, j, k) /= rhs; //scorro per la profontita'

	return a;
}

/**
 * Operator overloading = (assignment)
 *
 * Perform the assignment between this object and another
 *
 * @return a reference to the receiver object
 */
Tensor& Tensor::operator=(const Tensor& other)
{
	int n = other.r * other.c * other.d;

	this->r = other.r;
	this->c = other.c;
	this->d = other.d;

	if (this->data != other.data) {
		if (this->r < 0 && this->c < 0 && this->d < 0)
			throw unknown_operation();

		if (data)
			this->~Tensor();

		data = new float[n];

		for (int i = 0; i < n; i++) {
			this->data[i] = other.data[i];
		}
	}
	return (*this);
}

/**
 * Random Initialization
 *
 * Perform a random initialization of the tensor
 *
 * @param mean The mean
 * @param std  Standard deviation
 */
void Tensor::init_random(float mean, float std) {
	if (data) {

		std::default_random_engine generator;
		std::normal_distribution<float> distribution(mean, std);

		for (int i = 0; i < r; i++) {
			for (int j = 0; j < c; j++) {
				for (int k = 0; k < d; k++) {
					this->operator()(i, j, k) = distribution(generator);
				}
			}
		}

	}
	else {
		throw(tensor_not_initialized());
	}
}

/**
 * Constant Initialization
 *
 * Perform the initialization of the tensor to a value v
 *
 * @param r The number of rows
 * @param c The number of columns
 * @param d The depth
 * @param v The initialization value
 */
void Tensor::init(int r, int c, int d, float v) {
	//Controllo che v sia maggiore o uguale a 0
	if (v >= 0) {
		//Controllo che le dimensioni corrispondano
		if (this->r == r && this->c == c && this->d == d)
			for (int i = 0; i < this->r; i++)
				for (int j = 0; j < this->c; j++)
					for (int k = 0; k < this->d; k++)
						(*this)(i, j, k) = v;
		else
			throw dimension_mismatch();
	}
	else
		throw unknown_operation();
}

/**
 * Tensor Clamp
 *
 * Clamp the tensor such that the lower value becomes low and the higher one become high.
 *
 * @param low Lower value
 * @param high Higher value
 */
void Tensor::clamp(float low, float high) {
	//Controllo che il tensore sia inizializzato
	if (this->r != 0 && this->c != 0 && this->d != 0) {
		for (int i = 0; i < this->r; i++)
			for (int j = 0; j < this->c; j++)
				for (int k = 0; k < this->d; k++)
					if ((*this)(i, j, k) < low)
						(*this)(i, j, k) = low;
					else if ((*this)(i, j, k) > high)
						(*this)(i, j, k) = high;
	}
	else
		throw tensor_not_initialized();
}

/**
 * Tensor Rescaling
 *
 * Rescale the value of the tensor following this rule:
 *
 * newvalue(i,j,k) = ((data(i,j,k)-min(k))/(max(k)-min(k)))*new_max
 *
 * where max(k) and min(k) are the maximum and minimum value in the k-th channel.
 *
 * new_max is the new value for the maximum
 *
 * @param new_max New maximum vale
 */
void Tensor::rescale(float new_max) {
	//Controllo che il tensore sia inizializzato
	if (this->r != 0 && this->c != 0 && this->d != 0) {
		for (int k = 0; k < this->d; k++) {
			int min = this->getMin(k); //Minimo sulla profondita' k
			int max = this->getMax(k); //Massimo sulla profondita' k

			//Controllo che il minimo e il massimo siano diversi
			if (min == max)
				throw unknown_operation();

			for (int i = 0; i < this->r; i++)
				for (int j = 0; j < this->c; j++)
					(*this)(i, j, k) = (((*this)(i, j, k) - min) / (max - min)) * new_max;
		}

	}
	else
		throw tensor_not_initialized();
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
	if (r == 0 || c == 0 || d == 0) throw tensor_not_initialized{};//lancio errore di vettore non inizializzato 

	Tensor rit{ r + pad_h * 2, c + pad_w * 2, d, 0 }; //creo un nuovo oggetto di tipo tensore copiando quello passato

	for (int i = pad_h; i < rit.r - pad_h; i++)      //scorro per le righe
		for (int j = pad_w; j < rit.c - pad_w; j++) //scorro per le colonne
			for (int k = 0; k < rit.d; k++)      //scorro per la profondita' 
				rit(i, j, k) = (*this)(i - pad_h, j - pad_w, k); //sostituisco nel tensore di ritorno il valore del tensore su cui e' stata chiamata la funzione
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
	Tensor rit(row_end - row_start, col_end - col_start, depth_end - depth_start, 0); //creo un nuovo oggetto di tipo tensore copiando quello passato

	if (r == 0 || c == 0 || d == 0)               throw tensor_not_initialized{};//lancio errore di vettore non inizializzato 
	else if (rit.c > c || rit.r > r || rit.d > d) throw dimension_mismatch{};   //lancio errore di dimensioni non corrispondenti 

	for (size_t i = row_start; i < row_end; i++)               //scorro per le righe
		for (size_t j = col_start; j < col_end; j++)          //scorro per le colonne
			for (size_t k = depth_start; k < depth_end; k++) //scorro per la profondita' 
				rit(i - row_start, j - col_start, k - depth_start) = (*this)(i, j, k); //sostituisco nel tensore di ritorno il valore del tensore su cui e' stata chiamata la funzione
	return rit;
}

/**
* Concatenate
*
* The function concatenates two tensors along a give axis
*
* Example: this is of size 10x5x6 and rhs is of 25x5x6
*
* if concat on axis 0 (row) the result will be a new Tensor of size 35x5x6//colonne e prof =
*
* if concat on axis 1 (columns) the operation will fail because the number  //righe e prof=
* of rows are different (10 and 25).
*
* In order to perform the concatenation is mandatory that all the dimensions
* different from the axis should be equal, other wise throw concat_wrong_dimension().
*
* @param rhs The tensor to concatenate with
* @param axis The axis along which perform the concatenation
* @return a new Tensor containing the result of the concatenation
*/
Tensor Tensor::concat(const Tensor& rhs, int axis)const {
	int row{ this->r }, col{ this->c }, dep{ this->d };

	switch (axis) {
		//se asse=0
	case 0:
		if (col == rhs.c && dep == rhs.d)
			//calcolo righe del nuovo tensore
			row += rhs.r;
		else
			throw concat_wrong_dimension();
		break;

		//se asse=1
	case 1:
		if (row == rhs.r && dep == rhs.d)
			//calcolo colonne del nuovo tensore
			col += rhs.c;
		else
			throw concat_wrong_dimension();
		break;

		//se asse=2
	case 2:
		if (row == rhs.r && col == rhs.c)
			//calcolo profondità del nuovo tensore
			dep += rhs.d;
		else
			throw concat_wrong_dimension();
		break;

	default:
		throw unknown_operation();
		break;
	}

	//creo nuovo tensore
	Tensor t(row, col, dep);

	for (int i = 0; i < row; i++)
		for (int j = 0; j < col; j++)
			for (int k = 0; k < dep; k++) {
				//se siamo nel primo tensore
				if (i < this->r && j < this->c && k < this->d)
					t(i, j, k) = (*this)(i, j, k);
				//altrimenti
				else if (i >= this->r)
					t(i, j, k) = rhs(i - this->r, j, k);
				else if (j >= this->c)
					t(i, j, k) = rhs(i, j - this->c, k);
				else if (k >= this->d)
					t(i, j, k) = rhs(i, j, k - this->d);
			}

	return t;
}

/**
* Convolution
*
* This function performs the convolution of the Tensor with a filter.
*
* The filter f must have odd dimensions and same depth.
*
* Remeber to apply the padding before running the convolution
*
* @param f The filter
* @return a new Tensor containing the result of the convolution
*
*
* P = (F-1)/2 p=numero pixel da agg      F F è la dimensione del filtro ->R/C
*/
Tensor Tensor::convolve(const Tensor& f) const {
	//se dimensione scorretta
	if (this->d != f.d)
		throw dimension_mismatch();
	//se di dimensione pari
	if (f.r % 2 == 0 || f.c % 2 == 0)
		throw filter_odd_dimensions();

	//calcolo altezza e larghezza del padding	
	int hpad = (f.r - 1) / 2;
	int wpad = (f.c - 1) / 2;
	//creo tensore su cui lavorare con applicato il padding
	Tensor p = this->padding(hpad, wpad);
	Tensor ris(this->r, this->c, this->d);
	//per scorrere blocchi in orizzontale
	for (int i = 0; i < p.c - f.c + 1; i++) {
		//per scorrere blocchi in verticale
		for (int j = 0; j < p.r - f.r + 1; j++) {
			//per scorrere blocchi in profondità
			for (int k = 0; k < this->d; k++) {
				float somma = 0;
				//per scorrere le sotto-matrici
				for (int r = 0; r < f.r; r++) {
					for (int c = 0; c < f.c; c++) {
						somma += p(r + i, j + c, k) * f(r, c, k);
					}
				}
				//pongo all'interno del tensore ris i risultati della convoluzione
				ris(i, j, k) = somma;
			}
		}
	}
	return ris;
}

/* UTILITY */

/**
 * Rows
 *
 * @return the number of rows in the tensor
 */
int Tensor::rows()const {
	return this->r;
}

/**
 * Cols
 *
 * @return the number of columns in the tensor
 */
int Tensor::cols()const {
	return this->c;
}

/**
 * Depth
 *
 * @return the depth of the tensor
 */
int Tensor::depth()const {
	return this->d;
}

/**
 * Get minimum
 *
 * Compute the minimum value considering a particular index in the third dimension
 *
 * @return the minimum of data( , , k)
 */
float Tensor::getMin(int k)const {
	//Controllo che k sia in range [0, d]
	if (k >= 0 && k < this->d) {
		float min = (*this)(0, 0, 0); //Minimo del campo k

		for (int i = 0; i < r; i++)
			for (int j = 0; j < c; j++)
				if ((*this)(i, j, k) < min)
					min = (*this)(i, j, k);

		return min;
	}
	else
		throw index_out_of_bound();
}

/**
 * Get maximum
 *
 * Compute the maximum value considering a particular index in the third dimension
 *
 * @return the maximum of data( , , k)
 */
float Tensor::getMax(int k)const {
	//Controllo che k sia in range [0, d]
	if (k >= 0 && k < this->d) {
		float max = (*this)(0, 0, 0); //Massimo del campo k

		for (int i = 0; i < r; i++)
			for (int j = 0; j < c; j++)
				if ((*this)(i, j, k) > max)
					max = (*this)(i, j, k);

		return max;
	}
	else
		throw index_out_of_bound();
}

/**
 * showSize
 *
 * shows the dimensions of the tensor on the standard output.
 *
 * The format is the following:
 * rows" x "colums" x "depth
 *
 */
void Tensor::showSize()const {
	std::cout << this->r << "x" << this->c << "x" << this->d << std::endl;
}

/* IOSTREAM */

/**
 * Operator overloading <<
 *
 * Use the overaloading of << to show the content of the tensor.
 *
 * You are free to chose the output format, btw we suggest you to show the tensor by layer.
 *
 * [..., ..., 0]
 * [..., ..., 1]
 * ...
 * [..., ..., k]
 */
ostream& operator<< (ostream& stream, const Tensor& obj) {
	if (obj.r == 0 || obj.c == 0 || obj.d == 0) stream << "Il tensore non e' inizializzato";

	for (int k = 0; k < obj.d; k++)//scorro per la profondita'sore non e' inizzializzato";
	{
		stream << "LIVELLO " << k << endl;//aggiungo un identificatore che mi definisce a che profondità sono
		for (int j = 0; j < obj.c; j++)//scorro per la colonne
		{
			for (int i = 0; i < obj.r; i++)//scorro per la righe
				stream << obj(i, j, k) << "\t";     //aggiungo il valore allo stream e lo tabulo per seguire le colonne

			stream << endl;//vado a capo ad ogni colonna
		}
		stream << endl;//vado a capo ad ogni livello di profondità e salto una riga
	}

	return stream;
}

/**
 * Reading from file
 *
 * Load the content of a tensor from a textual file.
 *
 * The file should have this structure: the first three lines provide the dimensions while
 * the following lines contains the actual data by channel.
 *
 * For example, a tensor of size 4x3x2 will have the following structure:
 * 4
 * 3
 * 2
 * data(0,0,0)
 * data(0,1,0)
 * data(0,2,0)
 * data(1,0,0)
 * data(1,1,0)
 * .
 * .
 * .
 * data(3,1,1)
 * data(3,2,1)
 *
 * if the file is not reachable throw unable_to_read_file()
 *
 * @param filename the filename where the tensor is stored
 */
void Tensor::read_file(string filename) {
	ifstream file{ filename };

	//Controllo di aver aperto il file
	if (!file)
		throw unable_to_read_file();

	//Elimina il vettore data se esiste
	this->~Tensor();

	//Leggo dal file le dimensioni del tensore
	file >> this->r >> this->c >> this->d;

	//Creo il vettore dei dati
	data = new float[this->r * this->c * this->d];

	//Scrivo i dati del tensore nel file
	for (int k = 0; k < this->d; k++)
		for (int i = 0; i < this->r; i++)
			for (int j = 0; j < this->c; j++)
				file >> this->operator()(i, j, k);

	//Controllo se e' terminato il file
	if (!file.eof())
		throw dimension_mismatch();

	file.close();
}

/**
 * Write the tensor to a file
 *
 * Write the content of a tensor to a textual file.
 *
 * The file should have this structure: the first three lines provide the dimensions while
 * the following lines contains the actual data by channel.
 *
 * For example, a tensor of size 4x3x2 will have the following structure:
 * 4
 * 3
 * 2
 * data(0,0,0)
 * data(0,1,0)
 * data(0,2,0)
 * data(1,0,0)
 * data(1,1,0)
 * .
 * .
 * .
 * data(3,1,1)
 * data(3,2,1)
 *
 * @param filename the filename where the tensor should be stored
 */
void Tensor::write_file(string filename) {
	ofstream file{ filename, ios::trunc };

	//Controllo che il vettore sia inizializzato
	if (this->r == 0 || this->c == 0 || this->d == 0)
		throw tensor_not_initialized();

	//Scrivo nel file le dimensioni del tensore
	file << this->r << std::endl;
	file << this->c << std::endl;
	file << this->d << std::endl;

	//Scrivo nel file i dati del tensore
	for (int k = 0; k < this->d; k++)
		for (int i = 0; i < this->r; i++)
			for (int j = 0; j < this->c; j++) {
				file << this->operator()(i, j, k);
				if (i != this->r - 1 && j != this->c - 1 && k != this->d - 1)
					file << std::endl;
			}

	file.close();
}