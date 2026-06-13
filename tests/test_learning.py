"""Tests for the ILP learning subpackage.

These tests cover the task builder and base types.  The actual Aleph/Popper
backends require external dependencies (SWI-Prolog, Clingo) and are tested
separately with integration tests.
"""

import tempfile
from pathlib import Path

from pyaligner.learning.base import ILPTask, LearnedProgram
from pyaligner.learning.task_builder import AlignmentTaskBuilder


class TestILPTask:
    def test_empty_task(self):
        task = ILPTask()
        assert task.background == []
        assert task.positive == []
        assert task.negative == []
        assert task.bias == []

    def test_task_with_data(self):
        task = ILPTask(
            background=["parent(ann,bob)."],
            positive=["grandparent(ann,carol)."],
            negative=["grandparent(bob,ann)."],
            bias=[":- modeh(1, grandparent(+person, +person))."],
        )
        assert len(task.background) == 1
        assert len(task.positive) == 1
        assert len(task.negative) == 1


class TestLearnedProgram:
    def test_empty_program(self):
        prog = LearnedProgram()
        assert prog.clauses == []
        assert prog.program_text == ""

    def test_program_text(self):
        prog = LearnedProgram(
            clauses=["similar(X) :- seq1(X,S), length(S,L), L > 5.", "similar(X) :- sw_score(X,S), S > 10."],
            score=0.95,
        )
        assert "similar(X)" in prog.program_text
        assert prog.score == 0.95


class TestAlignmentTaskBuilder:
    def test_add_positive_pair(self):
        builder = AlignmentTaskBuilder()
        builder.add_positive_pair([1, 2, 3], [4, 5, 6], label="similar")
        task = builder.build()

        assert len(task.positive) == 1
        assert "similar(p1)." in task.positive[0]
        assert any("seq1(p1,[1,2,3])." in bg for bg in task.background)
        assert any("seq2(p1,[4,5,6])." in bg for bg in task.background)

    def test_add_negative_pair(self):
        builder = AlignmentTaskBuilder()
        builder.add_negative_pair([1, 2], [10, 11])
        task = builder.build()

        assert len(task.negative) == 1
        assert any("seq1" in bg for bg in task.background)

    def test_multiple_pairs(self):
        builder = AlignmentTaskBuilder()
        builder.add_positive_pair([1, 2], [3, 4])
        builder.add_positive_pair([5, 6], [7, 8])
        builder.add_negative_pair([1, 2], [10, 11])
        task = builder.build()

        assert len(task.positive) == 2
        assert len(task.negative) == 1
        # 2 background facts per pair (seq1 + seq2), 3 pairs = 6
        assert len(task.background) == 6

    def test_background_fact(self):
        builder = AlignmentTaskBuilder()
        builder.add_background_fact("amino_acid(1, a).")
        task = builder.build()
        assert "amino_acid(1, a)." in task.background

    def test_default_aleph_bias(self):
        builder = AlignmentTaskBuilder()
        builder.use_default_alignment_bias_aleph("similar")
        task = builder.build()
        assert any("modeh" in b for b in task.bias)
        assert any("determination" in b for b in task.bias)

    def test_default_popper_bias(self):
        builder = AlignmentTaskBuilder()
        builder.use_default_alignment_bias_popper("similar")
        task = builder.build()
        assert any("head_pred" in b for b in task.bias)
        assert any("body_pred" in b for b in task.bias)

    def test_write_files(self):
        builder = AlignmentTaskBuilder()
        builder.add_positive_pair([1, 2, 3], [4, 5, 6])
        builder.add_negative_pair([1, 2, 3], [10, 11, 12])
        builder.use_default_alignment_bias_popper()

        with tempfile.TemporaryDirectory() as tmpdir:
            task = builder.write_files(Path(tmpdir), name="test_task")

            # Check Aleph-format files.
            assert (Path(tmpdir) / "test_task.b").exists()
            assert (Path(tmpdir) / "test_task.f").exists()
            assert (Path(tmpdir) / "test_task.n").exists()

            # Check Popper-format files.
            assert (Path(tmpdir) / "bk.pl").exists()
            assert (Path(tmpdir) / "exs.pl").exists()
            assert (Path(tmpdir) / "bias.pl").exists()

            # Verify content.
            exs_content = (Path(tmpdir) / "exs.pl").read_text()
            assert "pos(similar(p1))" in exs_content
            assert "neg(similar(p2))" in exs_content

            bk_content = (Path(tmpdir) / "bk.pl").read_text()
            assert "seq1(p1,[1,2,3])." in bk_content

    def test_settings(self):
        builder = AlignmentTaskBuilder()
        builder.set_setting("i", "3")
        builder.set_setting("noise", "0")
        task = builder.build()
        assert task.settings == {"i": "3", "noise": "0"}

    def test_fluent_api(self):
        """Builder methods should return self for chaining."""
        builder = AlignmentTaskBuilder()
        result = (
            builder
            .add_positive_pair([1], [2])
            .add_negative_pair([3], [4])
            .add_background_fact("fact(a).")
            .add_bias(":- modeh(1, test(+x)).")
            .set_setting("i", "2")
        )
        assert result is builder


class TestAlephLearner:
    """Unit tests for AlephLearner that don't require SWI-Prolog."""

    def test_parse_output_from_rules(self):
        from pyaligner.learning.aleph import AlephLearner

        raw = """
[Rule 1]
similar(A) :-
    seq1(A,B), length(B,C), C > 5.

[Rule 2]
similar(A) :-
    sw_score(A,B), B > 10.

"""
        clauses = AlephLearner._parse_output(raw, Path("/nonexistent"))
        assert len(clauses) == 2
        assert "similar(A)" in clauses[0]

    def test_invalid_mode_raises(self):
        import pytest

        with pytest.raises(ValueError, match="Unknown induce_mode"):
            from pyaligner.learning.aleph import AlephLearner
            AlephLearner(induce_mode="invalid_mode")

    def test_valid_modes(self):
        from pyaligner.learning.aleph import AlephLearner

        for mode in AlephLearner.VALID_MODES:
            learner = AlephLearner(induce_mode=mode)
            assert learner.induce_mode == mode
