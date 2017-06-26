#pragma once

#include <Basis.h>
#include <Vector.h>
#include <Matrix.h>
#include <memory.h>

#include "Shape.h"

namespace GeometricalSolid{

	class Solid{
				
	public:
		Solid(std::unique_ptr<Shape> shape);
		~Solid();
		
		Solid(Solid&& other) = default;
		Solid(const Solid& other) = delete;

		Solid& operator=(Solid&& other) = default;
		Solid& operator=(const Solid& other) = delete;
		
		const Shape* Shape() const;
		const GeometricalSpaceObjects::Basis<double>& Basis() const;
		const GeometricalSpaceObjects::Vector<double>& Velocity() const;
		const GeometricalSpaceObjects::Vector<double>& AngularVelocity() const;
		const GeometricalSpaceObjects::Vector<double>& Force() const;
		const GeometricalSpaceObjects::Vector<double>& Momentum() const;
		GeometricalSpaceObjects::Matrix<double> Inertia() const;
		double Mass() const;

		void Shape(std::unique_ptr<GeometricalSolid::Shape> shape);
		void Basis(const GeometricalSpaceObjects::Basis<double> & basis);
		void Velocity(const GeometricalSpaceObjects::Vector<double> & v);
		void AngularVelocity(const GeometricalSpaceObjects::Vector<double>& w);
		void Force(const GeometricalSpaceObjects::Vector<double>& f);
		void Momentum(const GeometricalSpaceObjects::Vector<double>& m);
		
		void AddForce(const GeometricalSpaceObjects::Vector<double> & v);
		void AddMomentum(const GeometricalSpaceObjects::Vector<double> m);

		void LoadFromIstream(std::istream & in);
		
		void UpdateVelocities(double dt);
		void UpdatePosition(double dt);
		
		void ResetForceAndMomemtum();
		
		void LockTranslation(bool xAxis, bool yAxis, bool zAxis);
		void LockRotation(bool xAxis, bool yAxis, bool zAxis);

		bool IsXTranslationLocked();
		bool IsYTranslationLocked();
		bool IsZTranslationLocked();
		
		bool IsXRotationLocked();
		bool IsYRotationLocked();
		bool IsZRotationLocked();
		
	private:
		GeometricalSpaceObjects::Basis<double> basis;
		GeometricalSpaceObjects::Vector<double> velocity,angularVelocity,force,momentum;
		GeometricalSpaceObjects::Vector<double> localMomentum;
		GeometricalSpaceObjects::Vector<int> lockVelocity;
		GeometricalSpaceObjects::Vector<int> lockAngularVelocity;
		std::unique_ptr<GeometricalSolid::Shape> shape;
	};
}

inline std::ostream & operator << (std::ostream & out, const GeometricalSolid::Solid& a){
	out << a.Basis() << "\n" << a.Velocity() << "\n" << a.AngularVelocity() << "\n" << a.Force() << "\n" << a.Momentum();
	return out;
}

inline std::istream & operator >> (std::istream & in, GeometricalSolid::Solid & a){
	a.LoadFromIstream(in);
	return in;
}

