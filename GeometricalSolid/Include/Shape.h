#pragma once

#include <Matrix.h>

namespace GeometricalSolid{

	class Shape{
	public:
		enum class Nature{
			Container,
			Particle,
			Both,
			Unknown
		};

		enum class Form{
			Sphere,
			Rectangle,
			Disk,
			Cone,
			Elbow,
			Unknown
		};
		
		virtual ~Shape() {}
		
		const double& Mass() const {
			return this->mass;
		}
		
		const GeometricalSpaceObjects::Matrix<double>& InvertedIntertia() const {
			return this->invertedInertia;
		}

		Shape::Nature Nature() const { return this->nature; }
		Shape::Form Form() const { return this->form; }
		
	protected:
		enum Nature nature;
		enum Form form;
		double mass;
		GeometricalSpaceObjects::Matrix<double> invertedInertia;
	};

}
