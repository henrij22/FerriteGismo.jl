import os
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
EXEMPTIONS_FILE = REPO_ROOT / ".github/.filename_exemptions"


def load_exemptions():
    """
    Load exemptions from .filename_exemptions.
    Supports:
      - Exact paths
      - Glob patterns (e.g. docs/*.MD)
      - Comments (#) and empty lines
    """
    exemptions = []
    if EXEMPTIONS_FILE.exists():
        with EXEMPTIONS_FILE.open() as f:
            for line in f:
                line = line.strip()
                if line and not line.startswith("#"):
                    exemptions.append(line)
    return exemptions

def is_exempt(path: Path, exemptions):
    rel_path = path.relative_to(REPO_ROOT).as_posix()
    for pattern in exemptions:
        if path.match(pattern) or rel_path == pattern:
            return True
    return False


def has_uppercase_filename(path: Path) -> bool:
    return any(c.isupper() for c in path.name)


def main():
    exemptions = load_exemptions()
    violations = []

    for root, _, files in os.walk(REPO_ROOT):
        if ".git" in root:
            continue

        for filename in files:
            path = Path(root) / filename
            rel_path = path.relative_to(REPO_ROOT)

            if is_exempt(path, exemptions):
                continue

            if has_uppercase_filename(path):
                violations.append(rel_path.as_posix())

    if violations:
        print("❌ Found files with uppercase letters in filenames:")
        for v in violations:
            print(f"  - {v}")
        sys.exit(1)

    print("✅ All filenames are lowercase.")


if __name__ == "__main__":
    main()