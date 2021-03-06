#include <iostream>
#include <string>

#include "dais_exc.h"
#include "tensor.h"
#include "libbmp.h"
#include "DAISGram.h"

using namespace std;

DAISGram::DAISGram() {}

DAISGram::~DAISGram() {}

/**
 * Load a bitmap from file
 *
 * @param filename String containing the path of the file
 */
void DAISGram::load_image(string filename) {
	BmpImg img = BmpImg();

	img.read(filename.c_str());

	const int h = img.get_height();
	const int w = img.get_width();

	data = Tensor(h, w, 3, 0.0);

	for (int i = 0; i < img.get_height(); i++) {
		for (int j = 0; j < img.get_width(); j++) {
			data(i, j, 0) = (float)img.red_at(j, i);
			data(i, j, 1) = (float)img.green_at(j, i);
			data(i, j, 2) = (float)img.blue_at(j, i);
		}
	}
}


/**
 * Save a DAISGram object to a bitmap file.
 *
 * Data is clamped to 0,255 before saving it.
 *
 * @param filename String containing the path where to store the image.
 */
void DAISGram::save_image(string filename) {

	data.clamp(0, 255);

	BmpImg img = BmpImg(getCols(), getRows());

	img.init(getCols(), getRows());

	for (int i = 0; i < getRows(); i++) {
		for (int j = 0; j < getCols(); j++) {
			img.set_pixel(j, i, (unsigned char)data(i, j, 0), (unsigned char)data(i, j, 1), (unsigned char)data(i, j, 2));
		}
	}

	img.write(filename);

}

/**
 * Get rows
 *
 * @return returns the number of rows in the image
 */
int DAISGram::getRows() { return data.rows(); }

/**
 * Get columns
 *
 * @return returns the number of columns in the image
 */
int DAISGram::getCols() { return data.cols(); }

/**
 * Get depth
 *
 * @return returns the number of channels in the image
 */
int DAISGram::getDepth() { return data.depth(); }

/**
 * Brighten the image
 *
 * It sums the bright variable to all the values in the image.
 *
 * Before returning the image, the corresponding tensor should be clamped in [0,255]
 *
 * @param bright the amount of bright to add (if negative the image gets darker)
 * @return returns a new DAISGram containing the modified object
 */
DAISGram DAISGram::brighten(float bright)
{
	DAISGram rit;
	rit.data = data + bright;//carico data nel valore di ritorno e gli sommo la luminosit?? 

	rit.data.clamp(0, 255);//aggiusto i valori usciti dal range tra 0 e 255 

	return rit;
}

/**
 * Create a grayscale version of the object
 *
 * A grayscale image is produced by substituting each pixel with its average on all the channel
 *
 * @return returns a new DAISGram containing the modified object
 */
DAISGram DAISGram::grayscale()
{
    float s{ 0.0 };
    Tensor a{ data.rows(),data.cols(), data.depth(), 0 };//creo il tensore che andra' nel valore di ritorno 
    DAISGram rit;
    rit.data = a;//carico a nel valore di ritorno 

    for (int i = 0; i < data.rows(); i++)     //scorro per le righe
        for (int j = 0; j < data.cols(); j++)//scorro per le colonne
        {
            s = 0.0;
            for (int k = 0; k < data.depth(); k++) s += data(i, j, k);//scorro per la profondita' del tensore per calcolare la media

            s /= data.depth();//divido per il numero di valori sommati

            for (int k = 0; k < data.depth(); k++)  rit.data(i, j, k) = s;//scorro per la profondita' del tensore per inserire il valore appena calcolato
        }

    return rit;
}


/**
 * Create a Warhol effect on the image
 * 
 * This function returns a composition of 4 different images in which the:
 * - top left is the original image
 * - top right is the original image in which the Red and Green channel are swapped
 * - bottom left is the original image in which the Blue and Green channel are swapped
 * - bottom right is the original image in which the Red and Blue channel are swapped
 *  
 * The output image is twice the dimensions of the original one.
 * 
 * @return returns a new DAISGram containing the modified object
 */
DAISGram DAISGram::warhol(){
    Tensor redToGreen(this->data); // Tensore con canali rosso e verde invertiti
    Tensor blueToGreen(this->data); // Tensore con canali blu e verde invertiti
    Tensor redToBlue(this->data); // Tensore con canali rosso e blu invertiti
    float swap;

    for(int i = 0; i < data.rows(); i++)
        for(int j = 0; j < data.cols(); j++){
            swap = redToGreen(i,j, 0);
            redToGreen(i, j, 0) = redToGreen(i, j, 1);
            redToGreen(i, j, 1) = swap;

            swap = blueToGreen(i, j, 2);
            blueToGreen(i, j, 2) = blueToGreen(i, j, 1);
            blueToGreen(i, j, 1) = swap;

            swap = redToBlue(i, j, 0);
            redToBlue(i, j, 0) = redToBlue(i, j, 2);
            redToBlue(i, j, 2) = swap;
        }
    
    Tensor top = data.concat(redToGreen, 1);
    Tensor bottom = blueToGreen.concat(redToBlue, 1);

    DAISGram result;
    result.data = top.concat(bottom, 0);

    return result;
}

/**
 * Sharpen the image
 * 
 * This function makes the image sharper by convolving it with a sharp filter
 * 
 * filter[3][3]
 *    0  -1  0
 *    -1  5 -1
 *    0  -1  0
 *  
 * Before returning the image, the corresponding tensor should be clamped in [0,255]
 * 
 * @return returns a new DAISGram containing the modified object
 */
DAISGram DAISGram::sharpen(){
    DAISGram result;

	// Creo il filtro
    Tensor filter(3, 3, data.depth());
    for(int k = 0; k < filter.depth(); k++){
        filter(0, 1, k) = -1;
        filter(1, 0, k) = -1;
        filter(1, 1, k) = 5;
        filter(1, 2, k) = -1;
        filter(2, 1, k) = -1;
    }

    result.data = data.convolve(filter);
    result.data.clamp(0, 255);

    return result;
}


/**
 * Emboss the image
 * 
 * This function makes the image embossed (a light 3D effect) by convolving it with an
 * embossing filter
 * 
 * filter[3][3]
 *    -2 -1  0
 *    -1  1  1
 *     0  1  2
 * 
 * Before returning the image, the corresponding tensor should be clamped in [0,255]
 *  
 * @return returns a new DAISGram containing the modified object
 */
DAISGram DAISGram::emboss(){
    DAISGram result;

	// Creo il filtro
    Tensor filter(3, 3, data.depth());
    for(int k = 0; k < filter.depth(); k++){
        filter(0, 0, k) = -2;
        filter(0, 1, k) = -1;
        filter(1, 0, k) = -1;
        filter(1, 1, k) = 1;
        filter(1, 2, k) = 1;
        filter(2, 1, k) = 1;
        filter(2, 2, k) = 2;
    }

    result.data = data.convolve(filter);
    result.data.clamp(0, 255);

    return result;
}

/**
* Smooth the image
*
* This function remove the noise in an image using convolution and an average filter
* of size h*h:
*
* c = 1/(h*h)
*
* filter[3][3]
*    c c c
*    c c c
*    c c c
*
* @param h the size of the filter
* @return returns a new DAISGram containing the modified object
*/
DAISGram DAISGram::smooth(int h) {
	if (h == 0)
		throw unknown_operation();
	// calcolo c come  1/(h*h)
	float c = static_cast<float>(1.0 / (h * h));
	//creo il filtro
	Tensor filtro(h, h, data.depth(), c);
	DAISGram ris;
	//applico la convolve
	ris.data = this->data.convolve(filtro);
	return ris;
}
/**
* Edges of an image
*
* This function extract the edges of an image by using the convolution
* operator and the following filter
*
*
* filter[3][3]
* -1  -1  -1
* -1   8  -1
* -1  -1  -1
*
* Remeber to convert the image to grayscale before running the convolution.
*
* Before returning the image, the corresponding tensor should be clamped in [0,255]
*
* @return returns a new DAISGram containing the modified object
*/
DAISGram DAISGram::edge() {
	DAISGram ris;
	//applico il grayscale
	ris = this->grayscale();
	//creo un tensore=filtro 
	Tensor filtro(3, 3, data.depth(), -1);
	//modifico i valori all'interno del filtro
	filtro(1, 1, 0) = 8;
	filtro(1, 1, 1) = 8;
	filtro(1, 1, 2) = 8;
	//applico la convoluzione 
	ris.data = ris.data.convolve(filtro);
	//applico il clamp 
	ris.data.clamp(0, 255);
	return ris;
}

/**
* Blend with another image
*
* This function generate a new DAISGram which is the composition
* of the object and another DAISGram object
*
* The composition follows this convex combination:
* results = alpha*this + (1-alpha)*rhs
*
* rhs and this obejct MUST have the same dimensions.
*
* @param rhs The second image involved in the blending
* @param alpha The parameter of the convex combination
* @return returns a new DAISGram containing the blending of the two images.
*/
DAISGram DAISGram::blend(const DAISGram& rhs, float alpha) {
	DAISGram ris;
	//se le immagini hanno dimensioni scorrette
	if (this->data.cols() != rhs.data.cols() || this->data.rows() != rhs.data.rows() || this->data.depth() != rhs.data.depth())
		throw dimension_mismatch();
	//se alpha ha dimensioni scorrette
	if (alpha < 0 || alpha>1)
		throw unknown_operation();
	//results = alpha * this + (1 - alpha) * rhs
	Tensor tot = this->data * alpha + rhs.data * (1 - alpha);
	ris.data = tot;
	return ris;
}

/**
 * Green Screen
 * 
 * This function substitutes a pixel with the corresponding one in a background image 
 * if its colors are in the surrounding (+- threshold) of a given color (rgb).
 * 
 * (rgb - threshold) <= pixel <= (rgb + threshold)
 * 
 * 
 * @param bkg The second image used as background
 * @param rgb[] The color to substitute (rgb[0] = RED, rgb[1]=GREEN, rgb[2]=BLUE) 
 * @param threshold[] The threshold to add/remove for each color (threshold[0] = RED, threshold[1]=GREEN, threshold[2]=BLUE) 
 * @return returns a new DAISGram containing the result.
 */  
DAISGram DAISGram::greenscreen(DAISGram & bkg, int rgb[], float threshold[]){
	if(this->data.rows() == bkg.data.rows() && this->data.cols() == bkg.data.cols() && this->data.depth() == bkg.data.depth()){
		Tensor t(this->data);
		DAISGram result;
		result.data = t;

		for(int i = 0; i < this->data.rows(); i++)
			for(int j = 0; j < this->data.cols(); j++){
				int k = 0;
				bool flag = true; // Booleano per uscire dal ciclo se almeno un pixel ?? fuori dal range

				while (k < this->data.depth() && flag){
					if(data(i, j, k) < (rgb[k] - threshold[k]) || data(i, j , k) > (rgb[k] + threshold[k]))
						flag = false;
					k++;
				}

				if(flag){
					for(k = 0; k < result.data.depth(); k++)
						result.data(i, j, k) = bkg.data(i, j, k);
				}
			}
		return result;
	}
	else
		throw dimension_mismatch();
}

/**
 * Equalize
 * 
 * Stretch the distribution of colors of the image in order to use the full range of intesities.
 * 
 * See https://it.wikipedia.org/wiki/Equalizzazione_dell%27istogramma
 * 
 * @return returns a new DAISGram containing the equalized image.
 */  
DAISGram DAISGram::equalize(){
    DAISGram r;
    r.data = Tensor(this->data);
	
    for(int k = 0; k<this->data.depth(); k++){
		Tensor a{256, 1, 1};

		// Creo istogramma con i valori del tensore this->data
        for(int i = 0; i<this->data.rows(); i++)
        	for(int j = 0; j<this->data.cols(); j++)
				a(this->data(i, j, k),0,0)++;

		// Calcolo cdf(v)
        for(int i = 1; i<a.rows(); i++) 
        	a(i, 0, 0) += a(i-1, 0, 0);

		float m{0.0};
		int i{0};
			
		// Cerco il minimo cdf(v)
		while(i < 256 && m == 0){
        	if(a(i, 0, 0) > 0){
				m = a(i, 0, 0);
			} 
			i++;
		}

		// Calcolo il divisore della funzione di equalizzazione
		float divider = (r.getRows()*r.getCols()) - m;

        /*        cdf(v)-cdfMin
           round( -------------- *255)
        	      (r*c)-cdfMin
        */
        for(int i = 0; i<this->data.rows(); i++)
            for(int j = 0; j<this->data.cols(); j++)
				r.data(i, j, k) = round(((a(this->data(i, j, k), 0, 0)-m)/divider)*255.0);

    	}
		return r;
}

/**
 * Generate Random Image
 *
 * Generate a random image from nois
 *
 * @param h height of the image
 * @param w width of the image
 * @param d number of channels
 * @return returns a new DAISGram containing the generated image.
 */
void DAISGram::generate_random(int h, int w, int d) {
	data = Tensor(h, w, d, 0.0);
	data.init_random(128, 50);
	data.rescale(255);
}
