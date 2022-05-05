import numpy as np

class ShipCalculations:
    """
    Mathematical functions
    """
    def __init__(self):
        self.AT = 850 #[m]
        self.rho_air = 1.25 #[kg/m^3]
        self.rho_water = 1025 #[kg/m^3]
        self.CD = 0.82 #dimentionless
        self.breadth =  32.2 #[m]
        self.waterlineLenght = 234.526 #[m]
        self.g = 9.81

    def get_vec_magnitude(self, vector):
        return np.linalg.norm(vector)

    def angle_between_vectors(self, vec1, vec2):
        """
        Returns angle in radians
        """                
        unit_vector_1 = vec1 / np.linalg.norm(vec1)
        unit_vector_2 = vec2 / np.linalg.norm(vec2)
        dot_product = np.dot(unit_vector_1, unit_vector_2)
        angle = np.arccos(dot_product)
        return angle
    
    def wind_resistance(self, rho, drag_coefficient, cross_section, wind_speed): 
        return 1/2*rho*drag_coefficient*cross_section*wind_speed**2

    def STAwave_formulae(self, rho, g, significant_wave_height, breadth, length): 
        return 1/16*rho*g*significant_wave_height**2*breadth*np.sqrt(breadth/length)
