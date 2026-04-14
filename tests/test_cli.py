"""Tests for the codon-topo CLI."""

import json

from click.testing import CliRunner

from codon_topo.cli import main


runner = CliRunner()


class TestCLIHelp:
    def test_main_help(self):
        result = runner.invoke(main, ["--help"])
        assert result.exit_code == 0
        assert "Codon Geometry" in result.output

    def test_each_subcommand_help(self):
        for cmd in [
            "filtration",
            "disconnections",
            "coloring",
            "bit-bias",
            "trna",
            "kras",
            "claims",
            "condlogit",
            "topology-avoidance-k43",
            "codonsafe",
            "metric-sensitivity",
            "mis-analysis",
            "phylo-sensitivity",
            "all",
        ]:
            result = runner.invoke(main, [cmd, "--help"])
            assert result.exit_code == 0, f"{cmd} --help failed: {result.output}"


class TestFiltration:
    def test_default_table(self):
        result = runner.invoke(main, ["filtration"])
        assert result.exit_code == 0
        assert "Filtration" in result.output

    def test_json_output(self):
        result = runner.invoke(main, ["filtration", "--json"])
        assert result.exit_code == 0
        data = json.loads(result.output)
        assert "twofold_pass" in data
        assert data["twofold_pass"] == 9

    def test_all_tables(self):
        result = runner.invoke(main, ["filtration", "--all-tables", "--json"])
        assert result.exit_code == 0
        data = json.loads(result.output)
        assert len(data) >= 20


class TestDisconnections:
    def test_default(self):
        result = runner.invoke(main, ["disconnections"])
        assert result.exit_code == 0
        assert "Ser" in result.output

    def test_json(self):
        result = runner.invoke(main, ["disconnections", "--json"])
        assert result.exit_code == 0
        data = json.loads(result.output)
        assert any(e["aa"] == "Ser" for e in data)

    def test_all_tables(self):
        result = runner.invoke(main, ["disconnections", "--all-tables", "--json"])
        assert result.exit_code == 0
        data = json.loads(result.output)
        assert len(data) > 20


class TestColoring:
    def test_small_run(self):
        result = runner.invoke(main, ["coloring", "--n", "50", "--seed", "135325"])
        assert result.exit_code == 0
        assert "Observed score" in result.output

    def test_json(self):
        result = runner.invoke(
            main, ["coloring", "--n", "50", "--seed", "135325", "--json"]
        )
        assert result.exit_code == 0
        data = json.loads(result.output)
        assert "observed_score" in data
        assert "p_value_conservative" in data


class TestBitBias:
    def test_default(self):
        result = runner.invoke(main, ["bit-bias"])
        assert result.exit_code == 0
        assert "Weighted p-value" in result.output

    def test_json(self):
        result = runner.invoke(main, ["bit-bias", "--json"])
        assert result.exit_code == 0
        data = json.loads(result.output)
        assert "chi2_p_value_weighted" in data


class TestTRNA:
    def test_default(self):
        result = runner.invoke(main, ["trna"])
        assert result.exit_code == 0
        assert "Elevated tRNA" in result.output

    def test_json(self):
        result = runner.invoke(main, ["trna", "--json"])
        assert result.exit_code == 0
        data = json.loads(result.output)
        assert data["n_with_elevated_trna"] >= 7


class TestKRAS:
    def test_offline(self):
        result = runner.invoke(main, ["kras", "--offline"])
        assert result.exit_code == 0
        assert "OFFLINE" in result.output

    def test_offline_json(self):
        result = runner.invoke(main, ["kras", "--offline", "--json"])
        assert result.exit_code == 0
        data = json.loads(result.output)
        assert data["mode"] == "offline"
        assert "G12V" in data["fano_predictions"]


class TestClaims:
    def test_default(self):
        result = runner.invoke(main, ["claims"])
        assert result.exit_code == 0
        assert "SUPPORTED" in result.output
        assert "REJECTED" in result.output

    def test_json(self):
        result = runner.invoke(main, ["claims", "--json"])
        assert result.exit_code == 0
        data = json.loads(result.output)
        assert len(data) >= 10
        statuses = {c["status"] for c in data}
        assert "supported" in statuses
        assert "rejected" in statuses


class TestTopologyAvoidance:
    def test_json(self):
        result = runner.invoke(main, ["topology-avoidance", "--json"])
        assert result.exit_code == 0, f"Failed: {result.output}"
        data = json.loads(result.output)
        assert "rate_observed" in data
        assert "fisher_p" in data
        assert data["rate_observed"] < data["rate_possible"]


class TestRhoSweep:
    def test_json(self):
        result = runner.invoke(main, ["rho-sweep", "--n", "50", "--json"])
        assert result.exit_code == 0, f"Failed: {result.output}"
        data = json.loads(result.output)
        assert len(data["per_rho"]) == 5


class TestDecompose:
    def test_json(self):
        result = runner.invoke(main, ["decompose", "--json"])
        assert result.exit_code == 0
        data = json.loads(result.output)
        assert "by_nucleotide_position" in data
        assert "top_aa_pairs" in data


class TestPerTable:
    def test_json(self):
        result = runner.invoke(main, ["per-table", "--n", "50", "--json"])
        assert result.exit_code == 0, f"Failed: {result.output}"
        data = json.loads(result.output)
        assert data["n_tables"] >= 20


class TestAll:
    def test_run_all(self, tmp_path):
        result = runner.invoke(
            main, ["all", "--output-dir", str(tmp_path), "--n", "50"]
        )
        assert result.exit_code == 0, f"Failed: {result.output}"
        assert "Done" in result.output
        json_files = list(tmp_path.glob("*.json"))
        assert len(json_files) >= 8
        for jf in json_files:
            data = json.loads(jf.read_text())
            assert data is not None

    def test_manuscript_stats_json(self, tmp_path):
        result = runner.invoke(
            main, ["all", "--output-dir", str(tmp_path), "--n", "50"]
        )
        assert result.exit_code == 0
        stats_path = tmp_path / "manuscript_stats.json"
        assert stats_path.exists(), "manuscript_stats.json not generated"
        stats = json.loads(stats_path.read_text())
        assert stats["_version"] == "0.3.0"
        assert "metrics" in stats
        assert "coloring" in stats
        assert "topology_avoidance_q6" in stats
        assert "topology_avoidance_k43" in stats
        assert "condlogit" in stats
        assert "trna" in stats
