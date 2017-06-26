#pragma once

#include "Shape.h"
#include <Matrix.h>

namespace GeometricalSolid{
	
	class Retangle: public Shape{
	public:
		Retangle(const double lenght, const double width, const double thickness, const double density):lenght(lenght), width(width), thickness(thickness), density(density) {
			this->init();
		}
		
		double Density() const { return this->density; }
		double Volume() const { return this->volume; }
		double Lenght() const { return this->lenght; }
		double Width const { return this->width; }
		double Thickness const { return this-> thickness; }
		
	private:
		void init() {
			this->nature = Shape::Nature::Container;
			this->form = Shape::Form::Rectangle;

			this->volume = this->lenght*this->width*this->thickness;
			this->mass = this->volume*this->density;			
			this->inertia.Element(0, 0, this->mass/12*(this->width*this->width+this->thickness*this->thickness));
			this->inertia.Element(1, 1, this->mass/12*(this->lenght*this->lenght+this->thickness*this->thickness));
			this->inertia.Element(2, 2, this->mass/12*(this->lenght*this->lenght+this->width*this->width));
			this->invertedInertia = inertia.MatrixInverse();
		}
		
		double lenght{1};
		double width{1};
		double thickness{1};
		double density{2500};
		double volume{1};
		double mass{2500};
		GeometricalSpaceObjects::Matrix<double> inertia;
	};
}
