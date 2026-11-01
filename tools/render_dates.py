"""Record drawing dates and synchronize captions without editing image bytes."""
from pathlib import Path
from datetime import date, datetime
import argparse
import hashlib
import json
import re
from urllib.parse import quote

ROOT = Path(__file__).resolve().parents[1]
REGISTRY = ROOT / 'provenance/render-dates.json'
IMAGE = re.compile(r'^(!\[[^\]\n]*\]\(([^)\n]+)\))[ \t]*$', re.MULTILINE)
CAPTION = re.compile(r'\n\n<!-- drawing-date -->\n.*?\n<!-- /drawing-date -->', re.DOTALL)


def image_path(document, raw):
    if raw.startswith(('https://', 'http://', 'data:')):
        return None
    path = (document.parent / raw).resolve()
    try:
        relative = path.relative_to(ROOT).as_posix()
    except ValueError as exc:
        raise ValueError(f'Image is outside this repository: {raw}') from exc
    if not path.is_file():
        raise ValueError(f'Missing image: {relative}')
    return relative


def caption(row):
    date.fromisoformat(row['date'])
    if row['kind'] == 'first_recorded':
        evidence = row['evidence']
        link = f"{evidence['repository']}/blob/{evidence['commit']}/{quote(evidence['path'])}"
        text = f"*First recorded drawing: **{row['date']}** ([record]({link})).*"
    elif row['kind'] == 'rendered':
        text = f"*Drawn: **{row['date']}**.*"
    else:
        raise ValueError(f"Unknown date kind: {row['kind']}")
    return '\n\n<!-- drawing-date -->\n' + text + '\n<!-- /drawing-date -->'


def synchronize(ledger, check):
    missing, changed = [], []
    updates = []
    count = 0
    documents = set(ledger['caption_documents'])
    for pattern in ledger.get('caption_discovery', []):
        documents.update(p.relative_to(ROOT).as_posix() for p in ROOT.glob(pattern)
                         if IMAGE.search(p.read_text(encoding='utf-8')))
    documents.difference_update(ledger.get('frozen_document_exclusions', []))
    for name in sorted(documents):
        document = (ROOT / name).resolve()
        document.relative_to(ROOT)
        original = document.read_text(encoding='utf-8')
        bare = CAPTION.sub('', original)

        def replace(match):
            nonlocal count
            path = image_path(document, match[2])
            if path is None:
                return match[0]
            if path not in ledger['images']:
                missing.append(path)
                return match[0]
            count += 1
            return match[0] + caption(ledger['images'][path])

        updated = IMAGE.sub(replace, bare)
        if updated != original:
            changed.append(name)
        if updated != original:
            updates.append((document, updated))
    if missing:
        raise ValueError('Register drawing dates first: ' + ', '.join(sorted(set(missing))))
    if check and changed:
        raise ValueError('Run sync for: ' + ', '.join(changed))
    if not check:
        for document, updated in updates:
            document.write_text(updated, encoding='utf-8')
    print(json.dumps({'dated_image_occurrences': count, 'documents': len(documents), 'updated': [] if check else changed}))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='command', required=True)
    record = sub.add_parser('record', help='Run immediately after creating a new image.')
    record.add_argument('image', help='Repository-relative image path.')
    record.add_argument('--date', help='Earlier drawing date, YYYY-MM-DD; requires evidence.')
    record.add_argument('--evidence', help='Render log or other source for an explicitly supplied date.')
    sub.add_parser('sync', help='Update dated captions in the README, gallery and current research pages.')
    sub.add_parser('check', help='Fail if an image or caption lacks its registered date.')
    args = parser.parse_args()
    ledger = json.loads(REGISTRY.read_text(encoding='utf-8'))
    if args.command != 'record':
        synchronize(ledger, check=args.command == 'check')
        return
    path = (ROOT / args.image).resolve()
    key = path.relative_to(ROOT).as_posix()
    if not path.is_file() or path.suffix.lower() not in {'.png', '.gif', '.jpg', '.jpeg', '.svg', '.webp'}:
        parser.error('Supply an existing image inside the repository.')
    if key in ledger['images']:
        print(f"Preserved first drawing date {ledger['images'][key]['date']}: {key}")
        return
    if args.date and not args.evidence:
        parser.error('--date requires --evidence; do not guess an earlier drawing date.')
    now = datetime.now().astimezone()
    drawn = date.fromisoformat(args.date) if args.date else now.date()
    if drawn > now.date():
        parser.error('A drawing date cannot be in the future.')
    ledger['images'][key] = {
        'date': drawn.isoformat(),
        'kind': 'rendered',
        'recorded_at': now.isoformat(),
        'evidence': args.evidence or 'Recorded immediately after the first successful render.',
        'sha256_at_registration': hashlib.sha256(path.read_bytes()).hexdigest(),
    }
    ledger['images'] = dict(sorted(ledger['images'].items()))
    REGISTRY.write_text(json.dumps(ledger, indent=2, ensure_ascii=False) + '\n', encoding='utf-8')
    print(f"Recorded drawing date {drawn.isoformat()}: {key}")


if __name__ == '__main__':
    main()
