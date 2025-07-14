from cluster_affinity.cli import get_default_args
from cluster_affinity.matrix import get_matrix_args
from pathlib import Path
from argparse import ArgumentError

import pytest

class TestClusterScript:
    parser = get_default_args("cluster Affinity test","CA testing enviroment")

    @pytest.mark.cli
    def test_cluster_script_base_call(self):
        args = self.parser.parse_args("t1.tre t2.tre".split())
        assert args.t1 == Path("t1.tre")
        assert args.t2 == Path("t2.tre")
        assert args.unrooted == False
        assert args.cli == False
        assert args.color_only == False
        assert args.filetype == None

    @pytest.mark.cli
    def test_cluster_script_arguments_unrooted(self):
        args = self.parser.parse_args("t1.tre t2.tre --unrooted".split())
        assert args.unrooted == True

    @pytest.mark.cli
    def test_cluster_script_arguments_color_only(self):
        args = self.parser.parse_args("t1.tre t2.tre --color_only".split())
        assert args.color_only == True

    @pytest.mark.cli
    def test_cluster_script_arguments_cli(self):
        args = self.parser.parse_args("t1.tre t2.tre --cli".split())
        assert args.cli == True

    @pytest.mark.cli
    def tselfest_cluster_script_arguments_filetype(self):
        args = self.parser.parse_args("t1.tre t2.tre --filetype newick".split())
        assert args.filetype == "newick"

class TestMatrixScript:
    parser = get_matrix_args()

    @pytest.mark.cli
    def test_matrix_script_base_call(self):
        args = self.parser.parse_args("t1.tre outfile".split())
        assert args.t == [Path("t1.tre")]
        assert args.outfile == Path("outfile")
        assert args.title == ""
        assert args.cost == None
        assert args.filetype == None
        assert args.average == False
        assert args.csv_output == None
        assert args.unrooted == False

    def test_matrix_script_variable_arguments(self):
        args = self.parser.parse_args("t1.tre t2.tre t3.tre outfile".split())
        assert args.t == [Path("t1.tre"),Path("t2.tre"),Path("t3.tre")]
        assert args.outfile == Path("outfile")

    def test_matrix_script_title(self):
        args = self.parser.parse_args("t1.tre outfile --title Test".split())
        assert args.title == "Test"


    @pytest.mark.cli
    def test_matrix_script_cost_spec(self):
        args = self.parser.parse_args("t1.tre outfile --cost cluster_support".split())
        assert args.cost == "cluster_support"

    @pytest.mark.cli
    def test_matrix_script_filetype(self):
        args = self.parser.parse_args("t1.tre outfile --filetype nexus".split())
        assert args.filetype == "nexus"

    @pytest.mark.cli
    def test_matrix_script_average(self):
        args = self.parser.parse_args("t1.tre outfile --average".split())
        assert args.average == True

    @pytest.mark.cli
    def test_matrix_script_unrooted(self):
        args = self.parser.parse_args("t1.tre outfile --unrooted".split())
        assert args.unrooted == True

    @pytest.mark.cli
    def test_matrix_script_csv_output(self):
        args = self.parser.parse_args("t1.tre outfile --csv_output test.csv".split())
        assert args.csv_output == Path("test.csv")

    @pytest.mark.cli
    def test_matrix_script_autoscale(self):
        args = self.parser.parse_args("t1.tre outfile --autoscale".split())
        assert args.autoscale == True
