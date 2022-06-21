import numpy as np


class Beamsplitter:
    """This class defines a beam splitter

    The matrix describing the mode transformation is:

        e^{iphi}*cos(theta)     -sin(theta)
        e^{iphi}*sin(theta)     cos(theta)

    Args:
        mode1 (int): the index of the first mode (the first mode is mode 1)
        mode2 (int): the index of the second mode
        theta (float): the beam splitter angle
        phi (float): the beam splitter phase
    """

    def __init__(self, mode1, mode2, theta, phi):
        self.mode1 = mode1
        self.mode2 = mode2
        self.theta = theta
        self.phi = phi

    def __repr__(self):
        repr = "\n Beam splitter between modes {} and {}: \n Theta angle: {:.2f} \n Phase: {:.2f}".format(
            self.mode1,
            self.mode2,
            self.theta,
            self.phi
            )
        return repr


class Interferometer:
    """This class defines an interferometer. 

    An interferometer contains an ordered list of variable beam splitters,
    represented here by BS_list. For BS in BS_list, BS[0] and BS[1] correspond to the labels of the two modes
    being interfered (which start at 1). The beam splitters implement the optical transformation defined in equation 1 of:
    Clements, William R., et al. "Optimal design for universal multiport interferometers." Optica 3.12 (2016): 1460-1465.
    This transformation is parametrized by BS[2] (theta) which determines the beam splitter reflectivity, 
    and by BS[3] (phi). The interferometer also contains a list of output phases described by output_phases.
    """

    def __init__(self):
        self.BS_list = []
        self.output_phases = []

    def add_BS(self, BS):
        """Adds a beam splitter at the output of the current interferometer

        Args:
            BS (Beamsplitter): a Beamsplitter instance
        """
        self.BS_list.append(BS)

    def add_phase(self, mode, phase):    
        """Use this to manually add a phase shift to a selected mode at the output of the interferometer
        
        Args:
            mode (int): the mode index. The first mode is mode 1
            phase (float): the real-valued phase to add
        """
        while mode > np.size(self.output_phases):
            self.output_phases.append(0)
        self.output_phases[mode-1] = phase

    def count_modes(self):          
        """Calculate number of modes involved in the transformation. 

        This is required for calculate_transformation and draw

        Returns:
            the number of modes in the transformation
        """
        highest_index = max([max([BS.mode1, BS.mode2]) for BS in self.BS_list])
        return highest_index

    def calculate_transformation(self):       
        """Calculate unitary matrix describing the transformation implemented by the interferometer
    
        Returns:
            complex-valued 2D numpy array representing the interferometer
        """
        N = int(self.count_modes())
        U = np.eye(N, dtype=np.complex_)

        for BS in self.BS_list:
            T = np.eye(N, dtype=np.complex_)
            T[BS.mode1 - 1, BS.mode1 - 1] = np.exp(1j * BS.phi) * np.cos(BS.theta)
            T[BS.mode1 - 1, BS.mode2 - 1] = -np.sin(BS.theta)
            T[BS.mode2 - 1, BS.mode1 - 1] = np.exp(1j * BS.phi) * np.sin(BS.theta)
            T[BS.mode2 - 1, BS.mode2 - 1] = np.cos(BS.theta)
            U = np.matmul(T,U)

        while np.size(self.output_phases) < N:  # Autofill for users who don't want to bother with output phases
            self.output_phases.append(0)

        D = np.diag(np.exp([1j * phase for phase in self.output_phases]))
        U = np.matmul(D,U)
        return U

    def draw(self, show_plot=True):  
        """Use matplotlib to make a drawing of the interferometer

        Args:
            show_plot (bool): whether to show the generated plot
        """

        import matplotlib.pyplot as plt
        plt.figure()
        N = self.count_modes()
        mode_tracker = np.zeros(N)

        for ii in range(N):
            plt.plot((-1, 0), (ii, ii), lw=1, color="blue")

        for BS in self.BS_list:
            x = np.max([mode_tracker[BS.mode1 - 1], mode_tracker[BS.mode2 - 1]])
            plt.plot((x+0.3, x+1), (N - BS.mode1, N - BS.mode2), lw=1, color="blue")
            plt.plot((x, x+0.3), (N - BS.mode1, N - BS.mode1), lw=1, color="blue")
            plt.plot((x, x+0.3), (N - BS.mode2, N - BS.mode2), lw=1, color="blue")
            plt.plot((x+0.3, x+1), (N - BS.mode2, N - BS.mode1), lw=1, color="blue")
            plt.plot((x+0.4, x+0.9), (N - (BS.mode2 + BS.mode1)/2, N - (BS.mode2 + BS.mode1)/2), lw=1, color="blue")
            reflectivity = "{:2f}".format(np.cos(BS.theta)**2)
            plt.text(x+0.9, N + 0.05 - (BS.mode2 + BS.mode1)/2, reflectivity[0:3], color="green", fontsize=7)

            plt.plot((x+0.15, x+0.15), (N+0.3-(BS.mode2 + BS.mode1)/2., N+0.7-(BS.mode2 + BS.mode1)/2.), lw=1, color="blue")
            circle = plt.Circle((x+0.15, N+0.5-(BS.mode2 + BS.mode1)/2.), 0.1, fill=False)
            plt.gca().add_patch(circle)
            phase = "{:2f}".format(BS.phi)
            if BS.phi > 0:
                plt.text(x+0.2, N+0.7-(BS.mode2 + BS.mode1)/2., phase[0:3], color="red", fontsize=7)
            else:
                plt.text(x+0.2, N+0.7-(BS.mode2 + BS.mode1)/2., phase[0:4], color="red", fontsize=7)
            if x > mode_tracker[BS.mode1-1]:
                plt.plot((mode_tracker[BS.mode1-1], x), (N-BS.mode1, N-BS.mode1), lw=1, color="blue")
            if x > mode_tracker[BS.mode2-1]:
                plt.plot((mode_tracker[BS.mode2-1], x), (N-BS.mode2, N-BS.mode2), lw=1, color="blue")
            mode_tracker[BS.mode1-1] = x+1
            mode_tracker[BS.mode2-1] = x+1

        max_x = np.max(mode_tracker)
        for ii in range(N):
            plt.plot((mode_tracker[ii], max_x+1), (N-ii-1, N-ii-1), lw=1, color="blue")
            while np.size(self.output_phases) < N:  # Autofill for users who don't want to bother with output phases
                self.output_phases.append(0)
            if self.output_phases[ii] != 0:
                plt.plot((max_x+0.5, max_x+0.5), (N-ii-1.2, N-ii-0.8), lw=1, color="blue")
                circle = plt.Circle((max_x+0.5, N-ii-1), 0.1, fill=False)
                plt.gca().add_patch(circle)
                phase = str(self.output_phases[ii])
                if BS.phi > 0:
                    plt.text(max_x+0.6, N-ii-0.8, phase[0:3], color="red", fontsize=7)
                else:
                    plt.text(max_x+0.6, N-ii-0.8, phase[0:4], color="red", fontsize=7)


        plt.text(max_x/2, -0.7, "green: BS reflectivity", color="green", fontsize=10)
        plt.text(max_x/2, -1.4, "red: phase shift", color="red", fontsize=10)
        plt.text(-1, N-0.3, "Light in", fontsize=10)
        plt.text(max_x+0.5, N-0.3, "Light out", fontsize=10)
        plt.gca().axes.set_ylim([-1.8, N+0.2])
        plt.axis("off")
        if show_plot:
            plt.show()


def triangle_decomposition(U):
    """Returns a triangular mesh of beam splitters implementing matrix U

    This code implements the decomposition algorithm in:
    Reck, Michael, et al. "Experimental realization of any discrete unitary operator." 
    Physical review letters 73.1 (1994): 58.

    Args:
        U (np.ndarray): complex-valued 2D numpy array representing the interferometer

    Returns:
        an Interferometer instance
    """
    I = Interferometer()
    N = int(np.sqrt(U.size))
    for ii in range(N-1):
        for jj in range(N-1-ii):
            modes = [N - jj - 1, N - jj]
            theta = np.arctan(abs(U[ii, N - 1 - jj] / U[ii, N - 2 - jj]))
            phi = -np.angle(-U[ii, N - 1 - jj] / U[ii, N - 2 - jj])
            invT = np.eye(N, dtype=np.complex_)
            invT[modes[0]-1, modes[0]-1] = np.exp(-1j * phi) * np.cos(theta)
            invT[modes[0]-1, modes[1]-1] = np.exp(-1j * phi) * np.sin(theta)
            invT[modes[1]-1, modes[0]-1] = -np.sin(theta)
            invT[modes[1]-1, modes[1]-1] = np.cos(theta)
            U = np.matmul(U, invT)
            I.BS_list.append(Beamsplitter(modes[0], modes[1], theta, phi))
    phases = np.diag(U)
    I.output_phases = [np.angle(i) for i in phases]
    return I


def square_decomposition(U):
    """Returns a rectangular mesh of beam splitters implementing matrix U

    This code implements the decomposition algorithm in:
    Clements, William R., et al. "Optimal design for universal multiport interferometers." 
    Optica 3.12 (2016): 1460-1465.

    Returns:
        an Interferometer instance
    """
    I = Interferometer()
    N = int(np.sqrt(U.size))
    left_T = []
    for ii in range(N-1):
        if np.mod(ii, 2) == 0:
            for jj in range(ii+1):
                modes = [ii - jj + 1, ii + 2 - jj]
                theta = np.arctan(abs(U[N-1-jj, ii-jj] / U[N-1-jj, ii-jj+1]))
                phi = np.angle(U[N-1-jj, ii-jj] / U[N-1-jj, ii-jj+1])
                invT = np.eye(N, dtype=np.complex_)
                invT[modes[0]-1, modes[0]-1] = np.exp(-1j * phi) * np.cos(theta)
                invT[modes[0]-1, modes[1]-1] = np.exp(-1j * phi) * np.sin(theta)
                invT[modes[1]-1, modes[0]-1] = -np.sin(theta)
                invT[modes[1]-1, modes[1]-1] = np.cos(theta)
                U = np.matmul(U, invT)
                I.BS_list.append(Beamsplitter(modes[0], modes[1], theta, phi))
        else:
            for jj in range(ii+1):
                modes = [N+jj-ii-1, N+jj-ii]
                theta = np.arctan(abs(U[N+jj-ii-1, jj] / U[N+jj-ii-2, jj]))
                phi = np.angle(-U[N+jj-ii-1, jj] / U[N+jj-ii-2, jj])
                T = np.eye(N, dtype=np.complex_)
                T[modes[0]-1, modes[0]-1] = np.exp(1j * phi) * np.cos(theta)
                T[modes[0]-1, modes[1]-1] = -np.sin(theta)
                T[modes[1]-1, modes[0]-1] = np.exp(1j * phi) * np.sin(theta)
                T[modes[1]-1, modes[1]-1] = np.cos(theta)
                U = np.matmul(T, U)         
                left_T.append(Beamsplitter(modes[0], modes[1], theta, phi))

    for BS in np.flip(left_T, 0):
        modes = [int(BS.mode1), int(BS.mode2)]
        invT = np.eye(N, dtype=np.complex_)
        invT[modes[0]-1, modes[0]-1] = np.exp(-1j * BS.phi) * np.cos(BS.theta)
        invT[modes[0]-1, modes[1]-1] = np.exp(-1j * BS.phi) * np.sin(BS.theta)
        invT[modes[1]-1, modes[0]-1] = -np.sin(BS.theta)
        invT[modes[1]-1, modes[1]-1] = np.cos(BS.theta)
        U = np.matmul(invT, U)
        theta = np.arctan(abs(U[modes[1]-1, modes[0]-1]/U[modes[1]-1, modes[1]-1]))
        phi = np.angle(U[modes[1]-1, modes[0]-1] / U[modes[1]-1, modes[1]-1])
        invT[modes[0]-1, modes[0]-1] = np.exp(-1j * phi) * np.cos(theta)
        invT[modes[0]-1, modes[1]-1] = np.exp(-1j * phi) * np.sin(theta)
        invT[modes[1]-1, modes[0]-1] = -np.sin(theta)
        invT[modes[1]-1, modes[1]-1] = np.cos(theta)
        U = np.matmul(U, invT) 
        I.BS_list.append(Beamsplitter(modes[0], modes[1], theta, phi))
    phases = np.diag(U)
    I.output_phases = [np.angle(i) for i in phases]
    return I


def random_unitary(N):
    """Returns a random NxN unitary matrix

    This code is inspired by Matlab code written by Toby Cubitt:
    http://www.dr-qubit.org/matlab/randU.m

    Args:
        N (int): dimension of the NxN unitary matrix to generate

    Returns:
        complex-valued 2D numpy array representing the interferometer
    """
    X = np.zeros([N, N], dtype=np.complex_)
    for ii in range(N):
        for jj in range(N):
            X[ii, jj] = (np.random.normal() + 1j * np.random.normal()) / np.sqrt(2)

    q, r = np.linalg.qr(X)
    r = np.diag(np.divide(np.diag(r), abs(np.diag(r))))
    U = np.matmul(q, r)

    return U
