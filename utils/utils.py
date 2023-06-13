import numpy as np
from chemfiles import Frame


def custom_distance(positions, cell_lengths):
    '''
    Calculate the distance matrix in periodic boundary conditions.
    :param positions: np.array. set of coordinates
    :param cell_lengths: simulation box size
    :return: distance matrix
    '''

    cell_lengths = np.array(cell_lengths)
    distance = np.abs(positions[:, None, :] - positions)
    distance = np.where(distance > 0.5 * cell_lengths, distance - cell_lengths, distance)

    return np.sqrt((distance ** 2).sum(axis=-1))


def custom_distance_to_point(positions, point, cell_lengths):
    '''
    Calculate the distance matrix in periodic boundary conditions.
    :param positions: np.array. set of coordinates
    :param cell_lengths: simulation box size
    :param point:
    :return: distance matrix
    '''

    cell_lengths = np.array(cell_lengths)
    distance = np.abs(positions[:, None, :] - point)
    distance = np.where(distance > 0.5 * cell_lengths, distance - cell_lengths, distance)

    return np.sqrt((distance ** 2).sum(axis=-1))

def calc_dr_correlation(dr_1, dr_2):
    '''

    :param dr_1:
    :param dr_2:
    :return:
    '''
    correlation = np.sum(dr_1 * dr_2)/np.sqrt(np.sum(dr_1*dr_1) * np.sum(dr_2*dr_2))
    return correlation


def calc_list_correlation(corr_calcer_list):
    '''

    :param corr_calcer_list:
    :return:
    '''
    corr_frame = 0
    corr_calculated_times_frame = 0
    for i in range(len(corr_calcer_list)):
        for j in range(i + 1, len(corr_calcer_list)):
            if corr_calcer_list[i].ready_to_calculate and corr_calcer_list[j].ready_to_calculate:
                corr_frame += np.array(calc_dr_correlation(corr_calcer_list[i].dr, corr_calcer_list[j].dr))
                corr_calculated_times_frame += 1
    return corr_frame, corr_calculated_times_frame


def linear(x, a, b):
    return a*x+b

def constant_func(x, const):
    return np.full((len(x), ), const)


if __name__ == '__main__':
    a = np.array([1, 1, 0])
    b = np.array([1, 1, 0])
    print(calc_dr_correlation(a, b))
