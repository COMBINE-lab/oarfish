use logsumexp;
use std::ops;

#[derive(PartialEq, Clone, Copy)]
pub(crate) struct LogSpace {
    v: f64,
}

const _LOG_EPSILON: f64 = -36.04365338911715;
pub(crate) const MIN_LOG_P: LogSpace = LogSpace::new_from_ln(_LOG_EPSILON);

impl LogSpace {
    pub const fn new_from_ln(x: f64) -> Self {
        Self { v: x }
    }

    pub fn new_from_linear(x: f64) -> Self {
        Self { v: x.ln() }
    }

    pub fn get_linear(&self) -> f64 {
        self.v.exp()
    }

    #[allow(dead_code)]
    pub fn get_ln(&self) -> f64 {
        self.v
    }

    pub fn max(&self, other: Self) -> Self {
        Self::new_from_ln(self.v.max(other.v))
    }
}

impl ops::Add<LogSpace> for LogSpace {
    type Output = LogSpace;

    fn add(self, rhs: LogSpace) -> Self::Output {
        LogSpace {
            v: logsumexp::LogAddExp::ln_add_exp(&self.v, rhs.v),
        }
    }
}

impl ops::AddAssign<LogSpace> for LogSpace {
    fn add_assign(&mut self, rhs: Self) {
        self.v = logsumexp::LogAddExp::ln_add_exp(&self.v, rhs.v);
    }
}

impl ops::Mul<LogSpace> for LogSpace {
    type Output = LogSpace;

    // because we are in log space / is -
    #[allow(clippy::suspicious_arithmetic_impl)]
    fn mul(self, rhs: LogSpace) -> Self::Output {
        LogSpace { v: self.v + rhs.v }
    }
}

impl ops::MulAssign<LogSpace> for LogSpace {
    // because we are in log space / is -
    #[allow(clippy::suspicious_op_assign_impl)]
    fn mul_assign(&mut self, rhs: Self) {
        self.v += rhs.v;
    }
}

impl ops::Div<LogSpace> for LogSpace {
    type Output = LogSpace;
    // because we are in log space / is -
    #[allow(clippy::suspicious_arithmetic_impl)]
    fn div(self, rhs: LogSpace) -> Self::Output {
        LogSpace { v: self.v - rhs.v }
    }
}

impl ops::DivAssign<LogSpace> for LogSpace {
    // because we are in log space / is -
    #[allow(clippy::suspicious_op_assign_impl)]
    fn div_assign(&mut self, rhs: Self) {
        self.v -= rhs.v;
    }
}
