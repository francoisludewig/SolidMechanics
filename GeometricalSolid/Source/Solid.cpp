#include "../Include/Solid.h"
#include <Quaternion.h>
#include <iostream>
#include <iomanip>

using namespace GeometricalSolid;
using namespace GeometricalSpaceObjects;

Solid::Solid(std::unique_ptr<GeometricalSolid::Shape> shape):basis(),velocity(0,0,0),angularVelocity(0,0,0),force(0,0,0),momentum(0,0,0), lockVelocity(1,1,1), lockAngularVelocity(1,1,1), shape(std::move(shape)) {}

Solid::~Solid() {}

const Shape* Solid::Shape() const { return this->shape.get();}

const Basis<double>& Solid::Basis() const{return basis;}

const Vector<double>& Solid::Velocity() const{return velocity;}

const Vector<double>& Solid::AngularVelocity() const{return angularVelocity;}

const Vector<double>& Solid::Force() const{return force;}

const Vector<double>& Solid::Momentum() const{return momentum;}

void Solid::Basis(const GeometricalSpaceObjects::Basis<double> & basis){this->basis = basis;}

void Solid::Velocity(const Vector<double> & v) {this->velocity = v;}

void Solid::AngularVelocity(const Vector<double>& w) {this->angularVelocity = w;}

void Solid::Force(const Vector<double>& f) {this->force = f;}

void Solid::Momentum(const Vector<double>& m) {this->momentum = m;}

void Solid::Shape(std::unique_ptr<GeometricalSolid::Shape> shape) {
	this->shape = std::move(shape);
}

void Solid::AddForce(const GeometricalSpaceObjects::Vector<double> & force) {
	this->force += force;
}

void Solid::AddMomentum(const GeometricalSpaceObjects::Vector<double> momentum){
	this->momentum += momentum;
}

void Solid::UpdateVelocities(double dt){
  velocity += dt*force/(this->shape->Mass());
	localMomentum = momentum;
	basis.Local(localMomentum);
	angularVelocity += dt*(this->shape->InvertedIntertia()*localMomentum);
	
	// Lock
	this->velocity.ComponantX(this->velocity.ComponantX() * this->lockVelocity.ComponantX());
	this->velocity.ComponantY(this->velocity.ComponantY() * this->lockVelocity.ComponantY());
	this->velocity.ComponantZ(this->velocity.ComponantZ() * this->lockVelocity.ComponantZ());
	
	this->angularVelocity.ComponantX(this->angularVelocity.ComponantX() * this->lockAngularVelocity.ComponantX());
	this->angularVelocity.ComponantY(this->angularVelocity.ComponantY() * this->lockAngularVelocity.ComponantY());
	this->angularVelocity.ComponantZ(this->angularVelocity.ComponantZ() * this->lockAngularVelocity.ComponantZ());
}

void Solid::UpdatePosition(double dt){
	basis += dt*velocity;
	basis *= Quaternion<double>(angularVelocity*dt);
}

void Solid::ResetForceAndMomemtum(){
	force.SetComponants(0,0,0);
	momentum.SetComponants(0,0,0);
}

void Solid::LockTranslation(bool xAxis, bool yAxis, bool zAxis) {
	this->lockVelocity.SetComponants(xAxis ? 0 : 1, yAxis ? 0 : 1, zAxis ? 0 : 1);
}

void Solid::LockRotation(bool xAxis, bool yAxis, bool zAxis) {
	this->lockAngularVelocity.SetComponants(xAxis ? 0 : 1, yAxis ? 0 : 1, zAxis ? 0 : 1);
}

bool Solid::IsXTranslationLocked() { return this->lockVelocity.ComponantX() == 0; }
bool Solid::IsYTranslationLocked() { return this->lockVelocity.ComponantY() == 0; }
bool Solid::IsZTranslationLocked() { return this->lockVelocity.ComponantZ() == 0; }

bool Solid::IsXRotationLocked() { return this->lockAngularVelocity.ComponantX() == 0; }
bool Solid::IsYRotationLocked() { return this->lockAngularVelocity.ComponantY() == 0; }
bool Solid::IsZRotationLocked() { return this->lockAngularVelocity.ComponantZ() == 0; }

void Solid::LoadFromIstream(std::istream & in){
	in >> this->basis >> this->velocity >> this->angularVelocity >> this->force >> this->momentum;
}

