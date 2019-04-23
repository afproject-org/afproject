import os
import subprocess
#import tempfile uncessecarry?
import numpy as np
from ete3 import Tree



def tree_replace_ids(tree, id2name_dic):
    for id, name in id2name_dic.items():
        tree = tree.replace(id, name)
    return tree

def quartet_dist(tree_str1, tree_str2, temp_dir='./', token='abc'):

    file_path1 = os.path.join(temp_dir, token+'1.newick')
    file_path2 = os.path.join(temp_dir, token+'2.newick')
    oh = open(file_path1, 'w')
    oh.write(tree_str1)
    oh.close()
    oh = open(file_path2, 'w')
    oh.write(tree_str2)
    oh.close()

    cmd = [
        'quartet_dist',
        '-v',
        file_path1,
        file_path2,
    ]    
    sp = subprocess.Popen(cmd, shell=False,
          stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    output, error = sp.communicate()
    os.remove(file_path1)
    os.remove(file_path2)
    if not error:
        line = output.decode().strip()
        sl = line.split()
        n_leaves = int(sl[0])
        n_quartets = int(sl[1]) # n_leaves choose 4
        qd = int(sl[2])
        nqd = float(sl[3])
        return qd, nqd
    else:
        print(tree_str1)
        print(tree_str2)
        print(error)
        

def run_fneigbor(phylip_str, temp_dir='./', token='a'):
    phyfilename = os.path.join(temp_dir, token+'.phy')
    phyfile = open(phyfilename, 'w')
    phyfile.write(phylip_str)
    phyfile.close()
    outfilename = os.path.join(temp_dir, token+'.txt')
    outtree = os.path.join(temp_dir, token+'.newick')
    cmd = [
        'fneighbor', '-datafile', phyfilename,
        '-outfile', outfilename, 
        '-outtreefile', outtree,
        '-noprogress', '-notreeprint'
    ]
    sp = subprocess.Popen(cmd, shell=False,
          stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    output, error = sp.communicate()

    l = []
    fh = open(outtree)
    for line in fh:
        l.append(line.strip())
    fh.close()
    os.remove(phyfilename)
    os.remove(outfilename)
    os.remove(outtree)
    return "".join(l)



def ete_compare(usr_tree_str, ref_tree_str, outgroup_id=None):
    qt = Tree(usr_tree_str)
    rt = Tree(ref_tree_str)
    if outgroup_id:
        qt.set_outgroup(outgroup_id)
        rt.set_outgroup(outgroup_id)

    res = qt.compare(rt, unrooted=False if outgroup_id else True)
    rf = res["rf"]
    max_rf = res["max_rf"]
    nrf = res["norm_rf"]
    effective_tree_size = res["effective_tree_size"]
    ref_edges_in_source = res["ref_edges_in_source"]
    source_edges_in_ref = res["source_edges_in_ref"]
    source_subtrees = res["source_subtrees"]
    common_edges = res["common_edges"]
    treeko_dist = res["treeko_dist"]
    source_edges = res["source_edges"]
    ref_edges = res["ref_edges"]
    return qt, rt, nrf, rf, max_rf, source_edges_in_ref, ref_edges_in_source, treeko_dist
     


def run_ftreedist(usr_tree_str, ref_tree_str, temp_dir='./', token='abc', n=0):
    ftree_path = os.path.join(temp_dir, token+'.newick')
    ohtree = open(ftree_path, 'w')
    ohtree.write('{}\n{}'.format(usr_tree_str, ref_tree_str))
    ohtree.close()

    freedistout_path = ftree_path + '.out'
    cmd = [
        'ftreedist', ftree_path,
        '-outfile', freedistout_path,
        '-dtype', 's',
        '-style', 's',
        '-noprogress', '-noroot', 'N'
    ]    
    sp = subprocess.Popen(cmd, shell=False,
          stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    output, error = sp.communicate()
    fh = open(freedistout_path)
    line = fh.readline()
    fh.close()
    os.remove(ftree_path)
    os.remove(freedistout_path)
    if line:
        sl = line.strip().split()
        rf = int(sl[2])
        maxrf = 2 * (n - 3) if n else False
        nrf = rf/maxrf if maxrf else False
        return nrf, rf, maxrf

def validate_newick_treeids(tree, tree_ids):
    for id in tree_ids:
        if id not in tree:
            return False
    return tree


if __name__ == '__main__':
    import io
    import unittest


    usr_tree_str = '((A,C),(D,(B,E)));' 
    ref_tree_str = '(((A,D),C),(B,E));'


    class TestTsvValidation(unittest.TestCase):
 
      def setUp(self):
          pass

      def test_ete_compare_ok1(self):
          results = ete_compare(usr_tree_str, usr_tree_str)
          self.assertTrue(results[0])
          self.assertIsInstance(results[0], Tree)
          self.assertEqual(results[2], 0)
          self.assertEqual(results[3], 0)

      def test_ete_compare_ok2(self):
          results = ete_compare(usr_tree_str, ref_tree_str)
          self.assertTrue(results[0])
          self.assertIsInstance(results[0], Tree)
          self.assertEqual(results[3], 2)


      def test_quartet_distance_ok1(self):
          results = quartet_dist(usr_tree_str, usr_tree_str)
          self.assertEqual(results[0], 0)
          self.assertEqual(results[1], 0)

      def test_quartet_distance_ok2(self):
          results = quartet_dist(usr_tree_str, ref_tree_str)
          self.assertEqual(results[0], 2)
          self.assertEqual(results[1], 0.4)


    unittest.main()