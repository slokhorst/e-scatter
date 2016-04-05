/**
 * @file src/sample-viewer/camera.h
 * @author Thomas Verduin <T.Verduin@tudelft.nl>
 * @author Sebastiaan Lokhorst <S.R.Lokhorst@tudelft.nl>
 */

#ifndef SAMPLE_VIEWER__CAMERA__HEADER_INCLUDED
#define SAMPLE_VIEWER__CAMERA__HEADER_INCLUDED

#include <common/constant.hh>

class camera
{
private:
	control* m_control;
	cpl::vector3 m_position = cpl::vector3(-1,-1,+1);
	double horAngle=0.25*constant::pi, verAngle=-0.25*constant::pi;
	cpl::vector3 m_direction = cpl::vector3(+1,+1,-1);
	double m_speed = 1;
public:
	camera(control& control) {
		m_control = &control;
	}
	void set_speed(double v) {
		m_speed = v;
	}
	void update(double dt) {
		std::pair<double,double> mmove = m_control->mouse_movement();
		horAngle += 0.01f*mmove.first;
		verAngle -= 0.01f*mmove.second;
		verAngle = std::max(verAngle,-constant::pi/2+0.001f);
		verAngle = std::min(verAngle,+constant::pi/2-0.001f);

		m_direction = cpl::vector3(
			std::sin(horAngle)*std::cos(verAngle),
			std::cos(horAngle)*std::cos(verAngle),
			std::sin(verAngle)
		);

		cpl::vector3 movement_local = dt*m_control->requested_movement();
		cpl::vector3 movement_global = cpl::vector3(
			movement_local.y*std::sin(horAngle)*std::cos(verAngle) + movement_local.x*std::cos(horAngle),
			movement_local.y*std::cos(horAngle)*std::cos(verAngle) - movement_local.x*std::sin(horAngle),
			movement_local.y*std::sin(verAngle)
		);
		m_position += m_speed*movement_global;
	}
	cpl::vector3 position() const {
		return m_position;
	}
	cpl::vector3 direction() const {
		return m_direction;
	}
};

#endif
