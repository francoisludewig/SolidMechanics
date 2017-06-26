#pragma once

#include "Shape.h"
#include <Matrix.h>

namespace GeometricalSolid{
	
	class Disk: public Shape{
	public:
		RetDiskangle(const double radius, const double thickness, const double density):radius(raidus), thickness(thickness), density(density) {
			this->init();
		}
		
		double Density() const { return this->density; }
		double Volume() const { return this->volume; }
		double Radius() const { return this->radius; }
		double Thickness const { return this-> thickness; }
		
	private:
		void init() {
			this->nature = Shape::Nature::Container;
			this->form = Shape::Form::Disk;
			this->volume = this->radius*this->radius*M_PI*this->thickness;
			this->mass = this->volume*this->density;
			this->inertia.Element(0, 0, this->mass*(this->radius*this->radius/4. + this->thickness*this->thickness/12.));
			this->inertia.Element(1, 1, this->mass*(this->radius*this->radius/4. + this->thickness*this->thickness/12.));
			this->inertia.Element(2, 2, this->mass/2*(this->radius*this->radius));
			this->invertedInertia = inertia.MatrixInverse();
		}
		
		double radius{1};
		double thickness{1};
		double density{2500};
		double volume{M_PI};
		double mass{2500*M_PI};
		GeometricalSpaceObjects::Matrix<double> inertia;
	};
}
