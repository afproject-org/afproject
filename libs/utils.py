import math
import numpy as np

from ete3 import Tree


class Matrix():
    """Distance matrix

    Attributes:
        id_list (list): list of sequence identifiers
        data (ndarray): 2-D array of distance values between pairs of seqs
    """

    def __init__(self, id_list, data):
        """
        Example:
            >>> id_list = ['seq1', 'seq2', 'seq3']
            >>> data
            [[ 0.          0.3531587   0.35509333]
             [ 0.3531587   0.          0.295394  ]
             [ 0.35509333  0.295394    0.        ]]
            >>> matrix = Matrix(id_list, data)

        """
        self.id_list = id_list
        self.data = data

    def normalize(self):
        """Normalize distance values to 0-1 range."""
        self.data /= self.max()

    def __iter__(self):
        """Iterate over a distance matrix."""
        size = self.data.shape[0]
        for i, j in itertools.combinations(range(size), 2):
            yield i, j, self.id_list[i], self.id_list[j], self.data[i][j]

    def format(self, decimal_places=7):
        lines = ["   {0}".format(len(self.id_list))]
        for i, line in enumerate(self.data):
            seqid = self.id_list[i][:10]
            l = ['{0:.{1}f}'.format(line[i], decimal_places)
                 for i in range(0, len(line))]
            l.insert(0, '{0: <10}'.format(seqid))
            lines.append("\n" + " ".join(l))
        return "".join(lines)

    def min(self):
        """Return minimum distance value in matrix"""
        return np.amin(self.data)

    def max(self):
        """Return maximum distance value in matrix"""
        return np.amax(self.data)

    def is_zero(self):
        """Return True if matrix contains only zeros"""
        return not np.count_nonzero(self.data)


def read_tsv_dpair(handle):
    d = {}
    for line in handle:
        if line.strip():
            try:
                line = line.decode('utf-8')
            except:
                pass
            line = line.replace('.fasta', '')
            sl = line.split()
            id1 = sl[0]
            id2 = sl[1]
            val = float(sl[2])
            if math.isnan(val):
                raise ValueError('{}\t{} distance is NaN.'.format(id1, id2))
            #if val < 0:
            #    raise ValueError('{}\t{} has nagative value.'.format(id1, id2))
            pair = tuple(sorted([id1, id2]))
            d[pair] = val
    return d

def validate_tsv(handle):
    try:
        dpair = read_tsv_dpair(handle)
    except ValueError as e:
        return str(e), None
    except IndexError:
        e = "Number of columns is less than 3."
        return e, None
    except:
        e = "File is not in tsv format."
        return e, None
    return None, dpair


def read_clean_tsv_data_dpair(handle):
    d = {}
    for line in handle:
        sl = line.split()
        d[(sl[0], sl[1])] = float(sl[2])
    return d


def read_phy_matrix(handle):
    for i, line in enumerate(handle, start=-1):
        try:
            line = line.decode('utf-8')
        except:
            pass
        line = line.replace('.fasta', '')
        if i == -1:
            n = int(line.strip())
            ids = ['' for i in range(n)]
            data = np.zeros([n, n])
        else:
            sl = line.split()
            id = sl[0]
            ids[i] = id
            for j in range(0, len(sl)-1):
                val = float(sl[j+1])
                if math.isnan(val):
                    raise ValueError('{}\t{} distance is NaN.'.format(i, j))
                #if val < 0:
                #    raise ValueError('{}\t{} has negative value.'.format(i, j))
                data[i][j] = val
                data[j][i] = val
    return Matrix(id_list=ids, data=data)


def validate_phy(handle):
    try:
        matrix = read_phy_matrix(handle)
    except ValueError as e:
        return str(e), None
    except:
        e = "File is not in Phylip format."
        return e, None
    return None, matrix


def read_newick_tree(handle):
    tree = handle.read()
    try:
        tree = tree.decode('utf-8')
    except:
        pass
    return tree


def validate_newick(handle):
    tree = read_newick_tree(handle)
    tree = tree.replace('\n', '')
    tree = tree.replace(' ', '')
    tree = tree.replace('.fasta', '')
    if tree:
        if not tree.endswith(';'):
            e = 'Newick tree should end with semicolon (;)'
            return e, None
        if tree.count('(') != tree.count(')'):
            e = 'Open and close parantheses don\'t match up'
            return e, None
        try:
            t = Tree(tree)
        except:
            e = 'Malformatted Newick format'
            return e, None
    else:
        e = 'No tree in a file.'
        return e, None
    return None, tree


def validate_tree(tree_newick, tree_ids):
    tree = Tree(tree_newick)
    ids = set([tax.name for tax in tree])
    if len(ids) != len(tree_ids):
        e = 'The number of taxa in your tree does not correspond reference tree.'
        return e, None
    for id in tree_ids:
        if id not in ids:
            e = 'Tree does not have all seq ids (e.g., {})'.format(id)
            return e, None
    return None, tree_newick


def dpair_to_matrix(dpair, id_lst=None):
    if not id_lst:
        ids = set([])
        for p in dpair:
            for id in p:
                ids.add(id)
        id_lst = list(ids)
        id_lst.sort()
    n_ids = len(id_lst)
    data = np.zeros([n_ids, n_ids])
    for i, id1 in enumerate(id_lst):
        for j, id2 in enumerate(id_lst):
            if i != j:
                pair = tuple(sorted([id1, id2]))
                data[i][j] = dpair[pair]
    return Matrix(id_list=id_lst, data=data)


def matrix_to_dpair(matrix):
    dpair = {}
    n_ids = len(matrix.id_list)
    for i in range(0, n_ids):
        for j in range(i+1, n_ids):
            pair = tuple(sorted([matrix.id_list[i], matrix.id_list[j]]))
            dpair[pair] = matrix.data[i][j]
    return dpair


def validate_dpair(dpair, seqid_list):
    n = len(seqid_list)
    no_combinations = (n * (n-1))/2
    l = []
    if no_combinations != len(dpair):
        e = "Number of sequences is not equal to {}.".format(n)
        return e, None
    for pair in dpair:
        e = "Input file has extra id: {}"
        if pair[0] not in seqid_list:
            return e.format(pair[0]), None
        if pair[1] not in seqid_list:
            return e.format(pair[1]), None
        l.append('{}\t{}\t{}'.format(pair[0], pair[1], dpair[pair]))
    return None, "\n".join(l)


def validate_matrix(matrix, seqid_list):
    n = len(seqid_list)
    l = []
    if len(matrix.id_list) != n:
        e = "Number of sequences is not equal to {}.".format(n)
        return e, None
    if set(seqid_list).difference(set(matrix.id_list)):
        e = "Input file is missing some sequences."
        return e, None
    clean_phy = matrix.format()
    return None, clean_phy



if __name__ == '__main__':
    import io
    import unittest
    seqid_list = ['Alpha', 'Beta', 'Gamma', 'Delta', 'Epsilon']
    tsv_str = "\n".join([
      "Alpha    Beta    1.000",
      "Alpha    Gamma   2.000",
      "Alpha    Delta   3.000",
      "Alpha    Epsilon 3.000",
      "Beta Gamma   2.000",
      "Beta Delta   3.000",
      "Beta Epsilon    3.000",
      "Gamma    Delta   3.000",
      "Gamma    Epsilon 3.000",
      "Delta    Epsilon 1.000",
    ])
    tsv_str_ValueError1 = "\n".join([
      "Alpha    Beta    1.000",
      "Alpha    Gamma   2.000",
      "Alpha    Delta   STRING",
      "Alpha    Epsilon 3.000",
      "Beta Gamma   2.000",
      "Beta Delta   3.000",
      "Beta Epsilon    3.000",
      "Gamma    Delta   3.000",
      "Gamma    Epsilon 3.000",
      "Delta    Epsilon 1.000",
    ])
    tsv_str_ValueError2 = "\n".join([
      "Alpha    Beta    1.000",
      "Alpha    Gamma   2.000",
      "Alpha    Delta 3.000",
      "Alpha    Epsilon nan",
      "Beta Gamma   2.000",
      "Beta Delta   3.000",
      "Beta Epsilon    3.000",
      "Gamma    Delta   3.000",
      "Gamma    Epsilon 3.000",
      "Delta    Epsilon 1.000",
    ])
    tsv_str_IndexError = "\n".join([
      "Alpha    Beta    1.000",
      "Alpha    Gamma   2.000",
      "Alpha    Delta",
      "Alpha    Epsilon",
      "Beta Gamma   2.000",
      "Beta Delta   3.000",
      "Beta Epsilon    3.000",
      "Gamma    Delta   3.000",
      "Gamma    Epsilon 3.000",
      "Delta    Epsilon 1.000",
    ])

    phy_str = "\n".join([
      "   5",
      "Alpha      0.000 1.000 2.000 3.000 3.000",
      "Beta       1.000 0.000 2.000 3.000 3.000",
      "Gamma      2.000 2.000 0.000 3.000 3.000",
      "Delta      3.000 3.000 3.000 0.000 1.000",
      "Epsilon    3.000 3.000 3.000 1.000 0.000",
    ])
    phylt_str = "\n".join([
      "   5",
      "Alpha      ",
      "Beta       1.000",
      "Gamma      2.000 2.000",
      "Delta      3.000 3.000 3.000",
      "Epsilon    3.000 3.000 3.000 1.000",
    ])
    phy_str_ValueError1 = "\n".join([
      "   5",
      "Alpha      0.000 1.000 2.000 3.000 3.000",
      "Beta       1.000 0.000 2.000 3.000 3.000",
      "Gamma      2.000 STRING 0.000 3.000 3.000",
      "Delta      3.000 3.000 3.000 0.000 1.000",
      "Epsilon    3.000 3.000 3.000 1.000 0.000",
    ])
    phy_str_ValueError2_nan = "\n".join([
      "   5",
      "Alpha      0.000 1.000 2.000 3.000 3.000",
      "Beta       1.000 0.000 2.000 3.000 3.000",
      "Gamma      2.000 nan 0.000 3.000 3.000",
      "Delta      3.000 3.000 3.000 0.000 1.000",
      "Epsilon    3.000 3.000 3.000 1.000 0.000",
    ])
    newick_str = '(B,(A,C,E),D);'
    newick_str_invalid1 = '(B,(A,C,E),D)'
    newick_str_invalid2 = '(B,(A,C,E),D'

    class TestTsvValidation(unittest.TestCase):
 
      def setUp(self):
          pass

      def test_validate_tsv_ok(self):
          tsv_fh = io.StringIO(tsv_str)
          err, d = validate_tsv(tsv_fh)
          self.assertIsNone(err)
          self.assertIsNotNone(d)
          self.assertEqual(len(d), 10)
          self.assertEqual(d[('Beta', 'Gamma')], 2.0)

      def test_validate_tsv_ValueError1(self):
          tsv_fh = io.StringIO(tsv_str_ValueError1)
          err, d = validate_tsv(tsv_fh)
          self.assertIsNotNone(err)
          self.assertIsInstance(err, str)
          self.assertIsNone(d)

      def test_validate_tsv_ValueError2_nan(self):
          tsv_fh = io.StringIO(tsv_str_ValueError2)
          err, d = validate_tsv(tsv_fh)
          self.assertIsNotNone(err)
          self.assertIsInstance(err, str)
          self.assertIsNone(d)

      def test_validate_tsv_IndexError(self):
          tsv_fh = io.StringIO(tsv_str_IndexError)
          err, d = validate_tsv(tsv_fh)
          self.assertIsNotNone(err)
          self.assertIsInstance(err, str)
          self.assertIsNone(d)


    class TestPhyValidation(unittest.TestCase):
 
      def setUp(self):
          pass

      def test_validate_phy_ok(self):
          phy_fh = io.StringIO(phy_str)
          err, matrix = validate_phy(phy_fh)
          self.assertIsNone(err)
          self.assertIsNotNone(matrix)
          self.assertIsInstance(matrix, Matrix)

      def test_validate_phylt_ok(self):
          phy_fh = io.StringIO(phylt_str)
          err, matrix = validate_phy(phy_fh)
          self.assertIsNone(err)
          self.assertIsNotNone(matrix)
          self.assertIsInstance(matrix, Matrix)

      def test_validate_tsv_ValueError1(self):
          phy_fh = io.StringIO(phy_str_ValueError1)
          err, matrix = validate_phy(phy_fh)
          self.assertIsNotNone(err)
          self.assertIsInstance(err, str)
          self.assertIsNone(matrix)

      def test_validate_tsv_ValueError2_nan(self):
          phy_fh = io.StringIO(phy_str_ValueError2_nan)
          err, matrix = validate_phy(phy_fh)
          self.assertIsNotNone(err)
          self.assertIsInstance(err, str)
          self.assertIsNone(matrix)


    class TestTsvPhyConversion(unittest.TestCase):
 
      def setUp(self):
          fh = io.StringIO(phy_str)
          self.matrix_sym = read_phy_matrix(fh)
          fh = io.StringIO(phylt_str)
          self.matrix_lt = read_phy_matrix(fh)
          fh = io.StringIO(tsv_str)
          self.dpair = read_tsv_dpair(fh)

      def test_validate_matsym_matlt_equality(self):
          self.assertEqual(self.matrix_sym.id_list, self.matrix_lt.id_list)
          self.assertTrue(np.array_equal(self.matrix_sym.data, self.matrix_lt.data))
         
      def test_tsv2phy(self):
          dpair = matrix_to_dpair(self.matrix_sym)
          self.assertEqual(dpair, self.dpair)
          dpair = matrix_to_dpair(self.matrix_lt)
          self.assertEqual(dpair, self.dpair)

      def test_phy2tsv(self):
          matrix = dpair_to_matrix(self.dpair)
          self.assertEqual(len(matrix.id_list), len(self.matrix_sym.id_list))
          self.assertEqual(len(matrix.id_list), len(self.matrix_sym.id_list))
          self.assertEqual(set(matrix.id_list), set(self.matrix_sym.id_list))
          self.assertEqual(set(matrix.id_list), set(self.matrix_lt.id_list))
          self.assertEqual(len(matrix.data), len(self.matrix_sym.data))
          self.assertEqual(len(matrix.data), len(self.matrix_lt.data))
          self.assertEqual(matrix.data.shape, self.matrix_lt.data.shape)

      def test_phy2tsv_withIds(self):
          matrix = dpair_to_matrix(self.dpair, self.matrix_sym.id_list)
          self.assertEqual(matrix.id_list, self.matrix_sym.id_list)       
          self.assertEqual(matrix.id_list, self.matrix_lt.id_list)
          self.assertTrue(np.array_equal(matrix.data, self.matrix_sym.data))
          self.assertTrue(np.array_equal(matrix.data, self.matrix_lt.data))


    class TestDpairValidation(unittest.TestCase):
 
      def setUp(self):
          fh = io.StringIO(tsv_str)
          self.dpair = read_tsv_dpair(fh)

      def test_validate_dpair_ok(self):
          err, clean_data = validate_dpair(self.dpair, seqid_list)
          self.assertIsNone(err)
          self.assertTrue(clean_data)
          self.assertIsInstance(clean_data, str)

      def test_validate_dpair_missingPair(self):
          del self.dpair[('Alpha', 'Beta')]
          err, clean_data = validate_dpair(self.dpair, seqid_list)
          self.assertIsNotNone(err)
          self.assertIsNone(clean_data)

      def test_validate_dpair_extraPair(self):
          self.dpair[('Alpha', 'Zeta')] = 1.0
          err, clean_data = validate_dpair(self.dpair, seqid_list)
          self.assertIsNotNone(err)
          self.assertIsNone(clean_data)


    class TestMatrixValidation(unittest.TestCase):
 
      def setUp(self):
          fh = io.StringIO(phy_str)
          self.matrix_sym = read_phy_matrix(fh)
          fh = io.StringIO(phylt_str)
          self.matrix_lt = read_phy_matrix(fh)

      def test_validate_matrix_ok(self):
          err, clean_data = validate_matrix(self.matrix_sym, seqid_list)
          self.assertIsNone(err)
          self.assertTrue(clean_data)
          self.assertIsInstance(clean_data, str)
          err, clean_data = validate_matrix(self.matrix_sym, seqid_list)
          self.assertIsNone(err)
          self.assertTrue(clean_data)
          self.assertIsInstance(clean_data, str)

      def test_validate_matrix_invalidIds(self):
          seqid_list = ['seq{}'.format(i) for i in range(5)]
          err, clean_data = validate_matrix(self.matrix_sym, seqid_list)
          self.assertIsNotNone(err)
          self.assertIsInstance(err, str)
          self.assertIsNone(clean_data)

      def test_validate_matrix_missingSeq(self):
          err, clean_data = validate_matrix(self.matrix_sym, seqid_list + ['Zeta'])
          self.assertIsNotNone(err)
          self.assertIsInstance(err, str)
          self.assertIsNone(clean_data)

      def test_validate_matrix_extraSeq(self):
          err, clean_data = validate_matrix(self.matrix_sym, seqid_list[:-1])
          self.assertIsNotNone(err)
          self.assertIsInstance(err, str)
          self.assertIsNone(clean_data)

      def test_validate_matrix_to_phylip(self):
          self.assertEqual(self.matrix_sym.format(decimal_places=3), phy_str)


    class TestTreeValidation(unittest.TestCase):
 
      def setUp(self):
          fh = io.StringIO(newick_str)
          self.tree_newick = fh
          fh = io.StringIO(newick_str_invalid1)
          self.tree_newick_invalid1 = fh
          fh = io.StringIO(newick_str_invalid2)
          self.tree_newick_invalid2 = fh

      def test_validate_tree_ok(self):
          err, clean_data = validate_newick(self.tree_newick)
          self.assertIsNone(err)
          self.assertTrue(clean_data)


      def test_validate_newick_ok(self):
          err, clean_data = validate_newick(self.tree_newick)
          self.assertIsNone(err)
          self.assertTrue(clean_data)

      def test_validate_newick_ok(self):
          err, clean_data = validate_newick(self.tree_newick)
          self.assertIsNone(err)
          self.assertTrue(clean_data)


      def test_validate_tree_ok(self):
          err, clean_data = validate_newick(self.tree_newick)
          self.assertIsNone(err)
          self.assertTrue(clean_data)
          err, clean_data = validate_tree(clean_data, ['A', 'B', 'C', 'D', 'E'])
          self.assertIsNone(err)
          self.assertTrue(clean_data)

      def test_validate_newick_invalid1(self):
          err, clean_data = validate_newick(self.tree_newick_invalid1)
          self.assertIsNotNone(err)
          self.assertIsInstance(err, str)
          self.assertFalse(clean_data)

      def test_validate_newick_invalid2(self):
          err, clean_data = validate_newick(self.tree_newick_invalid2)
          self.assertIsNotNone(err)
          self.assertIsInstance(err, str)
          self.assertFalse(clean_data)

      def test_validate_tree_invalid(self):
          err, clean_data = validate_newick(self.tree_newick)
          self.assertIsNone(err)
          self.assertTrue(clean_data)
          err, clean_data = validate_tree(clean_data, ['A', 'B', 'C', 'D'])
          self.assertIsNotNone(err)
          self.assertIsInstance(err, str)
          self.assertFalse(clean_data)

    unittest.main()