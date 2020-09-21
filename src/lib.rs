use std::f64::consts::PI;

pub type Vector2 = nalgebra::Vector2<f64>;

const TAU: f64 = 2.0 * PI;

const RECUSION_LIMIT: u32 = 8;
const FLOAT_EPSILON: f64 = 1.19209290e-7;
const PATH_DISTANCE_EPSILON: f64 = 1.0;

const CURVE_ANGLE_TOERANCE_EPSILON: f64 = 0.01;
const M_ANGLE_TOLERANCE: f64 = 0.0;
const M_CUSP_LIMIT: f64 = 0.0;

#[inline]
fn clamp_angle(x: f64) -> f64 {
    if x >= PI {
        TAU - x
    } else {
        x
    }
}

pub fn adaptive_bezier_curve(
    start: Vector2,
    c1: Vector2,
    c2: Vector2,
    end: Vector2,
    scale: f64,
) -> Vec<Vector2> {
    let distance_tolerance = (PATH_DISTANCE_EPSILON / scale).powi(2);
    let mut sample_points = Vec::new();
    sample_points.push(start);
    adaptive_bezier_curve_impl(
        start,
        c1,
        c2,
        end,
        &mut sample_points,
        distance_tolerance,
        0,
    );
    sample_points.push(end);
    sample_points
}

pub fn adaptive_bezier_curve_impl(
    p1: Vector2,
    p2: Vector2,
    p3: Vector2,
    p4: Vector2,
    points: &mut Vec<Vector2>,
    distance_tolerance: f64,
    level: u32,
) {
    if level > RECUSION_LIMIT {
        return;
    }
    let p12 = (p1 + p2) / 2.0;
    let p23 = (p2 + p3) / 2.0;
    let p34 = (p3 + p4) / 2.0;
    let p123 = (p12 + p23) / 2.0;
    let p234 = (p23 + p34) / 2.0;
    let p1234 = (p123 + p234) / 2.0;
    if level > 0 {
        // Enforce subdivision first time
        // Try to approximate the full cubic curve by a single straight line
        let d = p4 - p1;
        let d2 = (p2 - p4).perp(&d).abs();
        let d3 = (p3 - p4).perp(&d).abs();
        if d2 > FLOAT_EPSILON && d3 > FLOAT_EPSILON {
            // Regular care
            if (d2 + d3).powi(2) <= distance_tolerance * d.norm_squared() {
                // If the curvature doesn't exceed the distanceTolerance value
                // we tend to finish subdivisions
                if M_ANGLE_TOLERANCE < CURVE_ANGLE_TOERANCE_EPSILON {
                    points.push(p1234);
                    return;
                }
                // Angle & cusp condition
                // TODO: replace here with angle
                let a23 = (p3.y - p2.y).atan2(p3.x - p2.x);
                let da1 = clamp_angle((a23 - (p2.y - p1.y).atan2(p2.x - p1.x)).abs());
                let da2 = clamp_angle(((p4.y - p3.y).atan2(p4.x - p3.x) - a23).abs());
                if da1 + da2 < M_ANGLE_TOLERANCE {
                    // Finally we can stop the recursion
                    points.push(p1234);
                    return;
                }
                if M_CUSP_LIMIT != 0.0 {
                    if da1 > M_CUSP_LIMIT {
                        points.push(p2);
                        return;
                    }
                    if da2 > M_CUSP_LIMIT {
                        points.push(p3);
                        return;
                    }
                }
            }
        } else if d2 > FLOAT_EPSILON {
            // P_1, P_3, P_4 are collinear, P_2 is considerable
            if d2 * d2 <= distance_tolerance * d.norm_squared() {
                if M_ANGLE_TOLERANCE < CURVE_ANGLE_TOERANCE_EPSILON {
                    points.push(p1234);
                    return;
                }
                // Angle condition
                let da1 = clamp_angle(
                    ((p3.y - p2.y).atan2(p3.x - p2.x) - (p2.y - p1.y).atan2(p2.x - p1.x)).abs(),
                );
                if da1 < M_ANGLE_TOLERANCE {
                    points.push(p2);
                    points.push(p3);
                    return;
                }
                if M_CUSP_LIMIT != 0.0 && da1 > M_CUSP_LIMIT {
                    points.push(p2);
                    return;
                }
            }
        } else if d3 > FLOAT_EPSILON {
            // P_1, P_2, P_4 are collinear, P_3 is considerable
            if d3 * d3 <= distance_tolerance * d.norm_squared() {
                if M_ANGLE_TOLERANCE < CURVE_ANGLE_TOERANCE_EPSILON {
                    points.push(p1234);
                    return;
                }
                // Angle condition
                let da1 = clamp_angle(
                    ((p4.y - p3.y).atan2(p4.x - p3.x) - (p3.y - p2.y).atan2(p3.x - p2.x)).abs(),
                );
                if da1 < M_ANGLE_TOLERANCE {
                    points.push(p2);
                    points.push(p3);
                    return;
                }
                if M_CUSP_LIMIT != 0.0 && da1 > M_CUSP_LIMIT {
                    points.push(p3);
                    return;
                }
            }
        } else {
            // Collinear case
            let d = p1234 - (p1 + p4) / 2.0;
            if d.norm_squared() < distance_tolerance {
                points.push(p1234);
                return;
            }
        }
    }
    // Continue subdivision
    adaptive_bezier_curve_impl(p1, p12, p123, p1234, points, distance_tolerance, level + 1);
    adaptive_bezier_curve_impl(p1234, p234, p34, p4, points, distance_tolerance, level + 1);
}

#[cfg(test)]
mod tests {
    use super::{adaptive_bezier_curve, Vector2, FLOAT_EPSILON};

    #[test]
    fn simple() {
        let answer = vec![
            Vector2::new(20.0, 20.0),
            Vector2::new(27.12921142578125, 32.740386962890625),
            Vector2::new(39.34417724609375, 56.408416748046875),
            Vector2::new(53.46435546875, 87.040771484375),
            Vector2::new(63.99200439453125, 111.28897094726562),
            Vector2::new(69.82025146484375, 123.60739135742188),
            Vector2::new(73.9102554321289, 130.69841384887695),
            Vector2::new(76.61632537841797, 134.36077499389648),
            Vector2::new(79.37641143798828, 137.14487075805664),
            Vector2::new(82.25093841552734, 139.02817916870117),
            Vector2::new(85.30033111572266, 139.98817825317383),
            Vector2::new(88.58501434326172, 140.00234603881836),
            Vector2::new(92.16541290283203, 139.04816055297852),
            Vector2::new(96.1019515991211, 137.10309982299805),
            Vector2::new(100.4550552368164, 134.1446418762207),
            Vector2::new(105.28514862060547, 130.15026473999023),
            Vector2::new(113.55682373046875, 122.16708374023438),
            Vector2::new(126.81915283203125, 107.68692016601562),
            Vector2::new(143.07708740234375, 88.65768432617188),
            Vector2::new(174.13818359375, 51.190185546875),
            Vector2::new(200.0, 20.0),
        ];
        let start = Vector2::new(20.0, 20.0);
        let c1 = Vector2::new(100.0, 159.0);
        let c2 = Vector2::new(50.0, 200.0);
        let end = Vector2::new(200.0, 20.0);
        let scale = 2.0;
        let output = adaptive_bezier_curve(start, c1, c2, end, scale);
        assert_eq!(output.len(), answer.len());
        for i in 0..answer.len() {
            assert_eq!((output[i].x - answer[i].x).abs() < FLOAT_EPSILON, true);
            assert_eq!((output[i].y - answer[i].y).abs() < FLOAT_EPSILON, true);
        }
    }
}
