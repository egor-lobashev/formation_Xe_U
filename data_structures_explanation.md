superdict:
{
    (bubble_size, ID_1):
        (
            [
                [t1, x1, y1, z1],
                [t2, x2, y2, z2],
                ...
                [tM, xM, yM, zM]        # M - last timestep of this bubble
            ],
            
            [ID_1, ID_2, ID_3, ... , ID_N]      # sorted. N - number of atoms in this bubble
        )
}


association_count:
{
    (size_1, size_2, timestep):  number_of_such_associations
}

C:
{
    (timestep, size):   number_of_such_bubbles
}

radii:
{
    size2: radius
}
(size1 is defined in the script association.py, n_min and n_max)