#!/usr/bin/env -S python -m unittest
import unittest
from collections.abc import Callable, Iterable
from functools import cache
from pathlib import Path

import molcube as mc
from molcube.pdbreader.pdbreader import PdbReaderProject
from molcube.pdbreader.glycosylation import Glycosylation, Node
from molcube.pdbreader.enums import GLYP_NODE_TYPE, PREDEFINED_GLYCAN

type RRString = str | mc.typed_dicts.ResidueRequest


@cache
def get_api_key() -> str:
    path = Path('api_access_key.txt')
    if not path.exists():
        raise FileNotFoundError(f"Create API key and save it in {path}")
    with path.open() as file:
        return file.read().strip()


def get_pdbreader() -> PdbReaderProject:
    molcube = mc.API('localhost', 8000)
    molcube.authenticate(api_token=get_api_key())

    return molcube.create_pdb_reader_project()


def delete_test_projects() -> None:
    molcube = mc.API('localhost', 8000)
    molcube.authenticate(api_token=get_api_key())
    test_names = ('test-mutation', 'test-phosphorylation-1kdx', 'test-phosphorylation-2klu',
                  'test-phosphorylation', 'test-protonation', 'test-ssbond',
                  'test-ssbond-default', 'test-staple', 'test-term',
                  'test-glycan-1', 'test-glycan-2', 'test-glycan-3', 'test-glycan-4')
    projects = molcube.search_projects(perPage=9999999)['projects']
    to_delete: list[str] = [project['pk'] for project in projects if project['title'] in test_names]
    molcube.delete_projects(to_delete)


#
# test glycans
#

def get_glycan_1() -> Glycosylation:
    # ============================================================
    # Example 1: Create a simple N-Glycan (2 sugars)
    # ============================================================
    # Target structure: bDGlcNAc(1->4)bDGlcNAc(1->)ASN
    #
    # ASN (root)
    #  └── bDGlcNAc (1->)
    #       └── bDGlcNAc (1->4)

    # Step 1: Create root node (PROTEIN type requires chain and resid)
    root = Node(
        name="ASN",                      # Residue name (must match actual PDB residue)
        type=GLYP_NODE_TYPE.PROTEIN,     # Root type for user-added glycosylation
        chain="PROT_A",                  # Chain index from PDB
        resid="160"                      # Residue ID from PDB
    )

    # Step 2: Create Glycosylation object
    gly = Glycosylation(root=root)  # chain_index will be auto-assigned

    print(f"Root node created: id={root.id}, name={root.name}, type={root.type}")
    print(f"Registry after root: {list(gly._registry.keys())}")

    # Step 3: Add first sugar (directly attached to root)
    sugar1 = Node(
        name="DGlcNAc",                  # Sugar name (case-insensitive)
        type=GLYP_NODE_TYPE.BETA         # Beta linkage (can also use "B")
    )
    # Note: resid is auto-assigned for sugar nodes (1, 2, 3, ...)

    sugar1_id = gly.add_sugar(
        node=sugar1,
        parent_id=0,      # 0 is the root node id
        site1="1",        # Child linkage site (usually "1" for sugars)
        site2=""          # Parent site (empty for root)
    )

    print(f"Sugar1 added: id={sugar1_id}, name={sugar1.name}, resid={sugar1.resid}")

    # Step 4: Add second sugar (attached to sugar1 at site 4)
    sugar2 = Node(
        name="DGlcNAc",
        type="B"  # String "B" also works for BETA
    )

    sugar2_id = gly.add_sugar(
        node=sugar2,
        parent_id=sugar1_id,  # Attach to sugar1
        site1="1",            # Child linkage site
        site2="4"             # Parent linkage site (1->4 linkage)
    )

    print(f"Sugar2 added: id={sugar2_id}, name={sugar2.name}, resid={sugar2.resid}")

    return gly


def get_glycan_2() -> Glycosylation:
    # ============================================================
    # Example 2: Branched Glycan Structure
    # ============================================================

    # Create root node
    root = Node(name="ASN", type=GLYP_NODE_TYPE.PROTEIN, chain="PROT_A", resid="315")
    gly_branched = Glycosylation(root=root)

    # Add first sugar (core GlcNAc)
    sugar1 = Node(name="DGlcNAc", type=GLYP_NODE_TYPE.BETA)
    sugar1_id = gly_branched.add_sugar(sugar1, parent_id=0, site1="1", site2="")

    # Add second sugar (core Man)
    sugar2 = Node(name="DMan", type=GLYP_NODE_TYPE.BETA)
    sugar2_id = gly_branched.add_sugar(sugar2, parent_id=sugar1_id, site1="1", site2="4")

    # Add branch 1: Man at site 6
    branch1 = Node(name="DMan", type=GLYP_NODE_TYPE.BETA)
    branch1_id = gly_branched.add_sugar(branch1, parent_id=sugar2_id, site1="1", site2="6")  # noqa: F841

    # Add branch 2: Man at site 3
    branch2 = Node(name="DMan", type=GLYP_NODE_TYPE.BETA)
    branch2_id = gly_branched.add_sugar(branch2, parent_id=sugar2_id, site1="1", site2="3")  # noqa: F841

    return gly_branched


def get_glycan_3() -> Glycosylation:
    # ============================================================
    # Example 3: Glycan with Modifications
    # ============================================================
    # this test is broken; the example code claims second set_modification()
    # should cause an exception, but no validation is done within
    # set_modification(); only add_glycosylation().

    root = Node(name="ASN", type=GLYP_NODE_TYPE.PROTEIN, chain="PROT_A", resid="160")
    gly_mod = Glycosylation(root=root)

    # Add sugar
    sugar = Node(name="DGlcNAc", type=GLYP_NODE_TYPE.BETA)
    sugar_id = gly_mod.add_sugar(sugar, parent_id=0, site1="1", site2="")

    # Add modification to sugar node
    # Parameters: node_id, site (key), modification value
    gly_mod.set_modification(node_id=sugar_id, key="6", value="S")  # Sulfate at site 6

    print("Glycan with modification:")

    # Note: This will fail - cannot add modification to root node
    try:
        gly_mod.set_modification(node_id=0, key="6", value="S")
        assert False, "Not supposed to get here"
    except ValueError:
        print("\n✗ Expected error: Root nodes cannot have modifications")

    return gly_mod


def get_glycan_4() -> Glycosylation:
    # ============================================================
    # Example 4: Using Predefined N-Glycan (M3)
    # ============================================================

    # Create root (ASN for N-glycans)
    root = Node(name="ASN", type=GLYP_NODE_TYPE.PROTEIN, chain="PROT_A", resid="160")
    gly_m3 = Glycosylation(root=root)

    # Apply M3 predefined glycan
    gly_m3.apply_predefined_glycan(PREDEFINED_GLYCAN.M3)  # or "M3"

    print("M3 N-glycan structure:")
    print(f"Total nodes: {len(gly_m3._registry)}")

    return gly_m3


class PdbReaderTest(unittest.TestCase):
    def test_mutation(self) -> None:
        pdbreader = get_pdbreader()
        self.assertTrue(pdbreader.create_project(title='test-mutation', ff='charmmff', customPdb='files/2klu.cif'))
        pdbreader.set_defaults()
        pdbreader.add_mutation(chain_id='PROT_A', resid='364', new_resname='ASN')  # GLY 364 -> ASN
        pdbreader.add_mutation(chain_id='PROT_A', resid='365', new_resname='ALA')  # PRO 365 -> ALA

        self.assertTrue(pdbreader.confirm_chains())
        self.assertTrue(pdbreader.model_pdb())
        pdbreader.download_pdb(f"{pdbreader.title}.pdb")

    def test_phosphorylation_1kdx(self) -> None:
        pdbreader = get_pdbreader()
        self.assertTrue(pdbreader.create_project(title='test-phosphorylation-1kdx', ff='charmmff', customPdb='files/1kdx.pdb'))
        pdbreader.set_defaults()
        pdbreader.add_phosphorylation(chain_id='PROT_B', resid='133', patch='SP2')

        self.assertTrue(pdbreader.confirm_chains())
        self.assertTrue(pdbreader.model_pdb())
        pdbreader.download_pdb(f"{pdbreader.title}.pdb")

    def test_phosphorylation_2klu(self) -> None:
        pdbreader = get_pdbreader()
        self.assertTrue(pdbreader.create_project(title='test-phosphorylation-2klu', ff='charmmff', customPdb='files/2klu.cif'))
        pdbreader.set_defaults()
        pdbreader.add_phosphorylation(chain_id='PROT_A', resid='394', patch='SP1')
        pdbreader.add_phosphorylation(chain_id='PROT_A', resid='415', patch='SP1')
        pdbreader.add_phosphorylation(chain_id='PROT_A', resid='431', patch='SP1')

        self.assertTrue(pdbreader.confirm_chains())
        self.assertTrue(pdbreader.model_pdb())
        pdbreader.download_pdb(f"{pdbreader.title}.pdb")

    def assert_ssbonds_exist(self, pdbreader: mc.pdbreader.pdbreader.PdbReaderProject,
                             ssbonds: Iterable[tuple[RRString, RRString]]) -> None:
        project_ssbonds = pdbreader._model_options.get('ssbond', [])
        for ssbond_formatted in ssbonds:
            residue1, residue2 = map(pdbreader._chain_res_unformat, ssbond_formatted)
            ssbond_unformatted = mc.typed_dicts.SsbondRequest(residue1=residue1, residue2=residue2)
            self.assertIn(member=ssbond_unformatted, container=project_ssbonds,
                          msg=f"Missing expected ssbond: {ssbond_unformatted}")

    def test_protonation(self) -> None:
        pdbreader = get_pdbreader()
        self.assertTrue(pdbreader.create_project(title='test-protonation', ff='charmmff', customPdb='files/6iyc.pdb'))
        pdbreader.set_defaults()

        # TODO: fix ligand and glycosylation settings so they don't need to be disabled
        pdbreader.toggle_chains_by_type(disable=['glycan', 'standaloneLigand'])
        self.assertTrue(pdbreader.confirm_chains())

        # PROT_A
        pdbreader.add_protonation(chain_id='PROT_A', resid='283', patch='ASPP')
        pdbreader.add_protonation(chain_id='PROT_A', resid='296', patch='GLUP')
        pdbreader.add_protonation(chain_id='PROT_A', resid='333', patch='GLUP')
        pdbreader.add_protonation(chain_id='PROT_A', resid='364', patch='GLUP')
        # PROT_B
        pdbreader.add_protonation(chain_id='PROT_B', resid='257', patch='ASPP')
        pdbreader.add_protonation(chain_id='PROT_B', resid='385', patch='ASPP')
        # PROT_C
        pdbreader.add_protonation(chain_id='PROT_C', resid='140', patch='ASPP')
        # disulfide bonds
        self.assert_ssbonds_exist(pdbreader, [
            ('PROT_A 50', 'PROT_A 62'),
            ('PROT_A 140', 'PROT_A 159'),
            ('PROT_A 230', 'PROT_A 248'),
            ('PROT_A 586', 'PROT_A 620'),
        ])

        pdbreader._model_options['glycosylation'].clear()

        self.assertTrue(pdbreader.model_pdb())
        pdbreader.download_pdb(f"{pdbreader.title}.pdb")

    def test_ssbond(self) -> None:
        pdbreader = get_pdbreader()
        self.assertTrue(pdbreader.create_project(title='test-ssbond', ff='charmmff', customPdb='files/4hg6.pdb'))
        pdbreader.set_defaults()

        # TODO: fix ligand and glycosylation settings so they don't need to be disabled
        pdbreader.toggle_chains_by_type(disable=['glycan', 'standaloneLigand'])

        self.assert_ssbonds_exist(pdbreader, [('PROT_B 163', 'PROT_B 430')])

        self.assertTrue(pdbreader.confirm_chains())
        self.assertTrue(pdbreader.model_pdb())
        pdbreader.download_pdb(f"{pdbreader.title}.pdb")

    def test_2hac_default_includes_ssbond(self) -> None:
        pdbreader = get_pdbreader()
        self.assertTrue(pdbreader.create_project(title='test-ssbond-default', ff='charmmff', customPdb='files/2hac.cif'))
        pdbreader.set_defaults()

        expected_ssbond = mc.typed_dicts.SsbondRequest(
            residue1={'chainIndex': 'PROT_A', 'resid': '2'},
            residue2={'chainIndex': 'PROT_B', 'resid': '2'},
        )

        ssbonds = pdbreader._model_options.get('ssbond', [])
        self.assertTrue(ssbonds, msg=f"Missing expected ssbond: {expected_ssbond}")
        self.assertTrue(len(ssbonds) == 1, msg=f"Too many ssbonds for 2hac: {ssbonds}")
        self.assertIn(member=expected_ssbond, container=ssbonds, msg="Expected "
                      f"{expected_ssbond} but got {ssbonds[0]}")

    def test_stapling(self) -> None:
        pdbreader = get_pdbreader()
        self.assertTrue(pdbreader.create_project(title='test-staple', ff='charmmff', customPdb='files/1ubq.pdb'))
        pdbreader.set_defaults()
        pdbreader.add_staple('RMETA3', 'PROT_A 1', 'PROT_A 3')

        self.assertTrue(pdbreader.confirm_chains())
        self.assertTrue(pdbreader.model_pdb())
        pdbreader.download_pdb(f"{pdbreader.title}.pdb")

    def test_terminal_patching(self) -> None:
        pdbreader = get_pdbreader()
        self.assertTrue(pdbreader.create_project(title='test-term', ff='charmmff', customPdb='files/1kdx.pdb'))
        pdbreader.set_defaults()
        pdbreader.set_terminal_patch(chain_id='PROT_A', nter='GLYP')  # default for cter
        pdbreader.set_terminal_patch(chain_id='PROT_B')  # default for both

        self.assertTrue(pdbreader.confirm_chains())
        self.assertTrue(pdbreader.model_pdb())
        pdbreader.download_pdb(f"{pdbreader.title}.pdb")

    def test_ligand(self) -> None:
        pdbreader = get_pdbreader()
        self.assertTrue(pdbreader.create_project(title='test-ligand', ff='charmmff', customPdb='files/2zff.cif'))
        pdbreader.set_defaults()

        assert pdbreader.validate_standalone_sdf(resname='53U')
        assert pdbreader.confirm_chains()

        pdbreader.get_nonstandard_sdf()
        assert pdbreader.model_pdb()

        pdbreader.download_pdb('myfile.pdb')

    def _glycan_test(self, getter: Callable[[], Glycosylation], title: str = '') -> None:
        pdbreader = get_pdbreader()
        self.assertTrue(pdbreader.create_project(title=title, ff='charmmff', customPdb='files/3pqr.cif'))
        pdbreader.set_defaults()

        glyc = getter()
        self.assertTrue(state := pdbreader.get_glycosylation_state(glyc))
        print(state)

        if getter is get_glycan_3:
            with self.assertRaises(ValueError):
                pdbreader.add_glycosylation(glyc)
        else:
            pdbreader.add_glycosylation(glyc)

        assert pdbreader.confirm_chains()
        assert pdbreader.model_pdb()

        pdbreader.download_project(f'{pdbreader.title}.pdb')

    def test_glycan_1(self) -> None:
        self._glycan_test(get_glycan_1, title='test-glycan-1')

    def test_glycan_2(self) -> None:
        self._glycan_test(get_glycan_2, title='test-glycan-2')

    def test_glycan_3(self) -> None:
        self._glycan_test(get_glycan_3, title='test-glycan-3')

    def test_glycan_4(self) -> None:
        self._glycan_test(get_glycan_4, title='test-glycan-4')

#
# TODO: other tests
#
#   PTMs: ???
#   ligands?:
#   heme?:
#   lipid tail (not implemented)
