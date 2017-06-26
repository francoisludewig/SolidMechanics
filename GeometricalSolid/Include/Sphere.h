#pragma once

#include "Shape.h"
#include <Matrix.h>

namespace GeometricalSolid{
	
	class Sphere: public Shape{
	public:
		Sphere():radius(1),density(1) {
			this->init();
		}
		
		Sphere(const double radius, const double density):radius(radius),density(density) {
			this->init();
		}
		
		~Sphere() {}
				
		double Radius() const { return this->radius; }
		double Density() const { return this->density; }
		double Volume() const { return this->volume; }
		
	private:
		void init() {
			this->nature = Shape::Nature::Particle;
			this->form = Shape::Form::Sphere;
			this->volume = 4./3.*M_PI*radius*radius*radius;
			this->mass = this->volume*this->density;
			this->inertia.Element(0, 0, 2./5.*this->mass*this->radius*this->radius);
			this->inertia.Element(1, 1, 2./5.*this->mass*this->radius*this->radius);
			this->inertia.Element(2, 2, 2./5.*this->mass*this->radius*this->radius);
			this->invertedInertia = inertia.MatrixInverse();
		}
		
		double radius;
		double density;
		double volume;
		GeometricalSpaceObjects::Matrix<double> inertia;
	};

}
