// @ts-check
import { defineConfig } from 'astro/config';
import starlight from '@astrojs/starlight';

// Published to the `gh-pages` branch and served at https://combine-lab.github.io/oarfish/
// (see .github/workflows/deploy-docs.yml). Links between pages include the `/oarfish` base.
export default defineConfig({
	site: 'https://combine-lab.github.io',
	base: '/oarfish',
	integrations: [
		starlight({
			title: 'oarfish',
			description:
				'Fast, accurate and versatile transcript quantification from long-read RNA-seq data.',
			favicon: '/favicon.svg',
			social: [
				{ icon: 'github', label: 'GitHub', href: 'https://github.com/COMBINE-lab/oarfish' },
			],
			editLink: {
				baseUrl: 'https://github.com/COMBINE-lab/oarfish/edit/main/docs/',
			},
			customCss: ['./src/styles/custom.css'],
			sidebar: [
				{
					label: 'Getting started',
					items: [
						{ label: 'Installation', slug: 'getting-started/installation' },
						{ label: 'Usage examples', slug: 'getting-started/usage-examples' },
					],
				},
				{
					label: 'Guides',
					items: [
						{ label: 'Quantification modes', slug: 'guides/quantification-modes' },
						{ label: 'Input', slug: 'guides/input' },
						{ label: 'Mapping backends', slug: 'guides/mapping-backends' },
						{ label: 'Coverage and truncation models', slug: 'guides/models' },
						{ label: 'Inferential replicates', slug: 'guides/inferential-replicates' },
						{ label: 'Single-cell mode', slug: 'guides/single-cell' },
					],
				},
				{
					label: 'Reference',
					items: [
						{ label: 'Command-line options', slug: 'reference/command-line' },
						{ label: 'Output files', slug: 'reference/output' },
						{
							label: 'Release notes',
							link: 'https://github.com/COMBINE-lab/oarfish/releases',
							attrs: { target: '_blank' },
						},
					],
				},
			],
		}),
	],
});
