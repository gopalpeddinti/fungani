import os
from types import SimpleNamespace

from fungani import cli, core
from fungani.cli import parse_args


def test_cli_has_main_entrypoint():
    assert hasattr(cli, "main")
    assert callable(cli.main)


def test_parser_default_outdir():
    parser = parse_args(["reference.fasta", "test.fasta"])
    assert parser.outdir is None


def test_parser_user_outdir():
    parser = parse_args(["reference.fasta", "test.fasta", "-o", "tmp"])
    assert parser.outdir is not None


def test_main_creates_missing_output_dir(tmp_path, monkeypatch):
    ref = tmp_path / "reference.fasta"
    test = tmp_path / "test.fasta"
    ref.write_text(">ref\nACGTACGT\n")
    test.write_text(">test\nACGTACGT\n")

    outdir = tmp_path / "results" / "fungani"
    monkeypatch.setattr(
        core,
        "make_blast_db",
        lambda pathname, filename: os.path.join(pathname, "db"),
    )
    monkeypatch.setattr(core, "run_async", lambda *args, **kwargs: [None])
    monkeypatch.setattr(core, "parse_results", lambda *args, **kwargs: ([], []))

    args = SimpleNamespace(
        outdir=str(outdir),
        mode="fwd",
        clean=False,
        size=16,
        overlap=2,
        cpus=1,
        threshold=80,
        percent=100,
        reference=str(ref),
        test=str(test),
        onepass=False,
    )

    result = core.main(args)

    assert result is None
    assert outdir.exists()
    assert (outdir / "fungani.log").exists()
