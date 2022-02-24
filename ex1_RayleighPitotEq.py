class RayleighPitotEquation():
    """ The Rayleigh Pitot Equation to calculate p_0,2/p_1 from M1.
    """    
    def __init__(self, gamma):
        self.gamma = gamma
        
    def __call__(self, M1):
        return ((self.gamma + 1)**2 * M1**2 / (4*self.gamma*M1**2 - 2*(self.gamma - 1)))**(self.gamma
                /(self.gamma-1)) * (1 - self.gamma + 2*self.gamma*M1**2) / (self.gamma + 1)

    def derivative(self, M1):
        """Calculate the first-order derivative at a given M1.

        Args:
            M1 (float).

        Returns:
            float: Value of the derivative.
        """        
        return (2*self.gamma*abs(self.gamma+1)**((2*self.gamma)/(self.gamma-1))*(2*M1**2-1)*abs(M1)**((2*self.gamma)
                /(self.gamma-1)))/((self.gamma+1)*M1*(4*self.gamma*M1**2-2*(self.gamma-1))**(self.gamma/(self.gamma-1)))
    
    def Nsolve(self, p_ratio, epsilon = 1e-5):
        """Use the Newton's method to approximate M1 at a given value of p_0,2/p_1

        Args:
            p_ratio (float): p_0,2/p_1.
            epsilon (float, optional): Maximum error. Defaults to 1e-5.
        
        Returns:
            M1_approx (float): The approximated value of M1 with an error <= epsilon.
        """        
                
        M1_approx = 1.5 # Initial guess of M1
        error = 1.0 # Randomly initialize the error
        i = 0 # Iteration
        
        while error > epsilon:
            i += 1
            M1_approx -= (self(M1_approx) - p_ratio)/self.derivative(M1_approx)
            error = abs(self(M1_approx) - p_ratio)
            print("Iteration {}: M1_{} = {}, error = {}".format(i, i, M1_approx, error))
        
        return M1_approx


if __name__ == '__main__':
    f = RayleighPitotEquation(1.4)
    print(f(3.5))
    print(f.derivative(3.5))
    print(f.Nsolve(f(3.5)))